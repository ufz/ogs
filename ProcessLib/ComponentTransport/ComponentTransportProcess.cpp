/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComponentTransportProcess.h"

#include <cassert>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/drop.hpp>

#include "BaseLib/RunTime.h"
#include "ChemistryLib/ChemicalSolverInterface.h"
#include "MathLib/LinAlg/Eigen/EigenTools.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/NumericalStability/FluxCorrectedTransport.h"
#include "NumLib/NumericalStability/NumericalStabilization.h"
#include "NumLib/ODESolver/NonlinearSystem.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/SurfaceFlux/SurfaceFlux.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/ComputeResiduum.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace ComponentTransport
{
ComponentTransportProcess::ComponentTransportProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ComponentTransportProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme,
    std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
        chemical_solver_interface,
    bool const is_linear,
    bool const ls_compute_only_upon_timestep_change)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux)),
      _chemical_solver_interface(std::move(chemical_solver_interface)),
      _asm_mat_cache{is_linear, use_monolithic_scheme},
      _ls_compute_only_upon_timestep_change{
          ls_compute_only_upon_timestep_change}
{
    if (ls_compute_only_upon_timestep_change)
    {
        // TODO move this feature to some common location for all processes.
        if (!is_linear)
        {
            OGS_FATAL(
                "Using the linear solver compute() method only upon timestep "
                "change only makes sense for linear model equations.");
        }

        WARN(
            "You specified that the ComponentTransport linear solver will do "
            "the compute() step only upon timestep change. This is an expert "
            "option. It is your responsibility to ensure that "
            "the conditions for the correct use of this feature are met! "
            "Otherwise OGS might compute garbage without being recognized. "
            "There is no safety net!");

        WARN(
            "You specified that the ComponentTransport linear solver will do "
            "the compute() step only upon timestep change. This option will "
            "only be used by the Picard non-linear solver. The Newton-Raphson "
            "non-linear solver will silently ignore this setting, i.e., it "
            "won't do any harm, there, but you won't observe the speedup you "
            "probably expect.");
    }

    auto residuum_name = [&](auto const& pv) -> std::string
    {
        std::string const& pv_name = pv.getName();
        if (pv_name == "pressure")
        {
            return "LiquidMassFlowRate";
        }
        if (pv_name == "temperature")
        {
            return "HeatFlowRate";
        }
        return pv_name + "FlowRate";
    };

    if (_use_monolithic_scheme)
    {
        int const process_id = 0;
        for (auto const& pv : _process_variables[process_id])
        {
            _residua.push_back(MeshLib::getOrCreateMeshProperty<double>(
                mesh, residuum_name(pv.get()), MeshLib::MeshItemType::Node, 1));
        }
    }
    else
    {
        for (auto const& pv : _process_variables)
        {
            _residua.push_back(MeshLib::getOrCreateMeshProperty<double>(
                mesh, residuum_name(pv[0].get()), MeshLib::MeshItemType::Node,
                1));
        }
    }
}

void ComponentTransportProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    int const mesh_space_dimension = _process_data.mesh_space_dimension;

    _process_data.mesh_prop_velocity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity",
        MeshLib::MeshItemType::Cell, mesh_space_dimension);

    _process_data.mesh_prop_porosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "porosity_avg",
        MeshLib::MeshItemType::Cell, 1);

    std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>
        transport_process_variables;
    if (_use_monolithic_scheme)
    {
        const int process_id = 0;
        for (auto pv_iter = std::next(_process_variables[process_id].begin());
             pv_iter != _process_variables[process_id].end();
             ++pv_iter)
        {
            transport_process_variables.push_back(*pv_iter);
        }
    }
    else
    {
        // All process variables but the pressure and optionally the temperature
        // are transport variables.
        for (auto const& pv :
             _process_variables |
                 ranges::views::drop(_process_data.isothermal ? 1 : 2))
        {
            transport_process_variables.push_back(pv[0]);
        }
    }

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh_space_dimension, mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data, transport_process_variables);

    if (_chemical_solver_interface && !_use_monolithic_scheme)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ComponentTransportLocalAssemblerInterface::setChemicalSystemID,
            _local_assemblers, _chemical_solver_interface->getElementIDs());

        _chemical_solver_interface->initialize();
    }

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh_space_dimension, getExtrapolator(), _local_assemblers,
            &ComponentTransportLocalAssemblerInterface::getIntPtDarcyVelocity));

    _secondary_variables.addSecondaryVariable(
        "molar_flux",
        makeExtrapolator(
            mesh_space_dimension, getExtrapolator(), _local_assemblers,
            &ComponentTransportLocalAssemblerInterface::getIntPtMolarFlux));
}

void ComponentTransportProcess::setInitialConditionsConcreteProcess(
    std::vector<GlobalVector*>& x, double const t, int const process_id)
{
    if (!_chemical_solver_interface)
    {
        return;
    }

    if (process_id != static_cast<int>(x.size() - 1))
    {
        return;
    }

    std::for_each(
        x.begin(), x.end(),
        [](auto const process_solution)
        { MathLib::LinAlg::setLocalAccessibleVector(*process_solution); });

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::initializeChemicalSystem,
        _local_assemblers, _chemical_solver_interface->getElementIDs(),
        dof_tables, x, t);
}

void ComponentTransportProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble ComponentTransportProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (_use_monolithic_scheme)
    {
        dof_tables.push_back(std::ref(*_local_to_global_index_map));
    }
    else
    {
        std::generate_n(
            std::back_inserter(dof_tables), _process_variables.size(),
            [&]() { return std::ref(*_local_to_global_index_map); });
    }

    _asm_mat_cache.assemble(t, dt, x, x_prev, process_id, M, K, b, dof_tables,
                            _global_assembler, _local_assemblers,
                            pv.getActiveElementIDs());

    // TODO (naumov) What about temperature? A test with FCT would reveal any
    // problems.
    if (process_id == _process_data.hydraulic_process_id)
    {
        return;
    }
    auto const matrix_specification = getMatrixSpecifications(process_id);

    NumLib::computeFluxCorrectedTransport(_process_data.stabilizer, t, dt, x,
                                          x_prev, process_id,
                                          matrix_specification, M, K, b);
}

void ComponentTransportProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ComponentTransportProcess.");

    if (_use_monolithic_scheme)
    {
        OGS_FATAL(
            "The AssembleWithJacobian() for ComponentTransportProcess is "
            "implemented only for staggered scheme.");
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    std::generate_n(std::back_inserter(dof_tables), _process_variables.size(),
                    [&]() { return std::ref(*_local_to_global_index_map); });

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_prev, process_id, M, K, b, Jac);

    // b is the negated residumm used in the Newton's method.
    // Here negating b is to recover the primitive residuum.
    transformVariableFromGlobalVector(b, 0, *_local_to_global_index_map,
                                      *_residua[process_id],
                                      std::negate<double>());
}

Eigen::Vector3d ComponentTransportProcess::getFlux(
    std::size_t const element_id,
    MathLib::Point3d const& p,
    double const t,
    std::vector<GlobalVector*> const& x) const
{
    std::vector<GlobalIndexType> indices_cache;
    auto const r_c_indices = NumLib::getRowColumnIndices(
        element_id, *_local_to_global_index_map, indices_cache);

    std::vector<std::vector<GlobalIndexType>> indices_of_all_coupled_processes{
        x.size(), r_c_indices.rows};
    auto const local_xs =
        getCoupledLocalSolutions(x, indices_of_all_coupled_processes);

    return _local_assemblers[element_id]->getFlux(p, t, local_xs);
}

void ComponentTransportProcess::solveReactionEquation(
    std::vector<GlobalVector*>& x, std::vector<GlobalVector*> const& x_prev,
    double const t, double const dt, NumLib::EquationSystem& ode_sys,
    int const process_id)
{
    // todo (renchao): move chemical calculation to elsewhere.
    if (_process_data.lookup_table && process_id == 0)
    {
        INFO("Update process solutions via the look-up table approach");
        _process_data.lookup_table->lookup(x, x_prev, _mesh.getNumberOfNodes());

        return;
    }

    if (!_chemical_solver_interface)
    {
        return;
    }

    // Sequential non-iterative approach applied here to split the reactive
    // transport process into the transport stage followed by the reaction
    // stage.
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    if (process_id == 0)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ComponentTransportLocalAssemblerInterface::setChemicalSystem,
            _local_assemblers, _chemical_solver_interface->getElementIDs(),
            dof_tables, x, t, dt);

        BaseLib::RunTime time_phreeqc;
        time_phreeqc.start();

        _chemical_solver_interface->setAqueousSolutionsPrevFromDumpFile();

        _chemical_solver_interface->executeSpeciationCalculation(dt);

        INFO("[time] Phreeqc took {:g} s.", time_phreeqc.elapsed());

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ComponentTransportLocalAssemblerInterface::
                postSpeciationCalculation,
            _local_assemblers, _chemical_solver_interface->getElementIDs(), t,
            dt);

        return;
    }

    auto const matrix_specification = getMatrixSpecifications(process_id);

    std::size_t matrix_id = 0u;
    auto& M = NumLib::GlobalMatrixProvider::provider.getMatrix(
        matrix_specification, matrix_id);
    auto& K = NumLib::GlobalMatrixProvider::provider.getMatrix(
        matrix_specification, matrix_id);
    auto& b =
        NumLib::GlobalVectorProvider::provider.getVector(matrix_specification);

    M.setZero();
    K.setZero();
    b.setZero();

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::assembleReactionEquation,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t, dt, M, K,
        b, process_id);

    // todo (renchao): incorporate Neumann boundary condition
    MathLib::finalizeMatrixAssembly(M);
    MathLib::finalizeMatrixAssembly(K);
    MathLib::finalizeVectorAssembly(b);

    auto& A = NumLib::GlobalMatrixProvider::provider.getMatrix(
        matrix_specification, matrix_id);
    auto& rhs =
        NumLib::GlobalVectorProvider::provider.getVector(matrix_specification);

    A.setZero();
    rhs.setZero();

    MathLib::finalizeMatrixAssembly(A);
    MathLib::finalizeVectorAssembly(rhs);

    // compute A
    MathLib::LinAlg::copy(M, A);
    MathLib::LinAlg::aypx(A, 1.0 / dt, K);

    // compute rhs
    MathLib::LinAlg::matMult(M, *x[process_id], rhs);
    MathLib::LinAlg::aypx(rhs, 1.0 / dt, b);

    using Tag = NumLib::NonlinearSolverTag;
    using EqSys = NumLib::NonlinearSystem<Tag::Picard>;
    auto& equation_system = static_cast<EqSys&>(ode_sys);
    equation_system.applyKnownSolutionsPicard(A, rhs, *x[process_id]);

    auto& linear_solver =
        _process_data.chemical_solver_interface->linear_solver;
    MathLib::LinAlg::setLocalAccessibleVector(*x[process_id]);
    linear_solver.solve(A, rhs, *x[process_id]);

    NumLib::GlobalMatrixProvider::provider.releaseMatrix(M);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(K);
    NumLib::GlobalVectorProvider::provider.releaseVector(b);
    NumLib::GlobalMatrixProvider::provider.releaseMatrix(A);
    NumLib::GlobalVectorProvider::provider.releaseVector(rhs);
}

void ComponentTransportProcess::computeSecondaryVariableConcrete(
    double const t,
    double const dt,
    std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev,
    int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_prev, process_id);

    if (!_chemical_solver_interface)
    {
        return;
    }

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::
            computeReactionRelatedSecondaryVariable,
        _local_assemblers, _chemical_solver_interface->getElementIDs());
}

void ComponentTransportProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev,
    const double t,
    const double dt,
    int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::postTimestep,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, x_prev, t,
        dt, process_id);

    if (!_surfaceflux)  // computing the surfaceflux is optional
    {
        return;
    }
    _surfaceflux->integrate(x, t, *this, process_id, _integration_order, _mesh,
                            pv.getActiveElementIDs());
}

void ComponentTransportProcess::preOutputConcreteProcess(
    const double t,
    double const dt,
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev,
    int const process_id)
{
    auto const matrix_specification = getMatrixSpecifications(process_id);

    auto M = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(
        matrix_specification);
    auto K = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(
        matrix_specification);
    auto b = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        matrix_specification);

    M->setZero();
    K->setZero();
    b->setZero();

    assembleConcreteProcess(t, dt, x, x_prev, process_id, *M, *K, *b);

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (_use_monolithic_scheme)
    {
        dof_tables.push_back(std::ref(*_local_to_global_index_map));
    }
    else
    {
        std::generate_n(
            std::back_inserter(dof_tables), _process_variables.size(),
            [&]() { return std::ref(*_local_to_global_index_map); });
    }

    BaseLib::RunTime time_residuum;
    time_residuum.start();

    if (_use_monolithic_scheme)
    {
        auto const residuum =
            computeResiduum(dt, *x[0], *x_prev[0], *M, *K, *b);
        for (std::size_t variable_id = 0; variable_id < _residua.size();
             ++variable_id)
        {
            transformVariableFromGlobalVector(
                residuum, variable_id, dof_tables[0], *_residua[variable_id],
                std::negate<double>());
        }
    }
    else
    {
        auto const residuum = computeResiduum(dt, *x[process_id],
                                              *x_prev[process_id], *M, *K, *b);
        transformVariableFromGlobalVector(residuum, 0, dof_tables[process_id],
                                          *_residua[process_id],
                                          std::negate<double>());
    }

    INFO("[time] Computing residuum flow rates took {:g} s",
         time_residuum.elapsed());
}
}  // namespace ComponentTransport
}  // namespace ProcessLib
