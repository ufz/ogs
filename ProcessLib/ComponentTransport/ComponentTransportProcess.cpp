/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComponentTransportProcess.h"

#include <cassert>

#include "BaseLib/RunTime.h"
#include "ChemistryLib/ChemicalSolverInterface.h"
#include "MathLib/LinAlg/Eigen/EigenTools.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/SurfaceFlux/SurfaceFlux.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
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
        chemical_solver_interface)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data)),
      _surfaceflux(std::move(surfaceflux)),
      _chemical_solver_interface(std::move(chemical_solver_interface))
{
}

void ComponentTransportProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    _process_data.mesh_prop_velocity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity",
        MeshLib::MeshItemType::Cell, mesh.getDimension());

    _process_data.mesh_prop_porosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "porosity_avg",
        MeshLib::MeshItemType::Cell, 1);

    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>
        transport_process_variables;
    if (_use_monolithic_scheme)
    {
        for (auto pv_iter = std::next(_process_variables[process_id].begin());
             pv_iter != _process_variables[process_id].end();
             ++pv_iter)
        {
            transport_process_variables.push_back(*pv_iter);
        }
    }
    else
    {
        for (auto pv_iter = std::next(_process_variables.begin());
             pv_iter != _process_variables.end();
             ++pv_iter)
        {
            transport_process_variables.push_back((*pv_iter)[0]);
        }
    }

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data,
        transport_process_variables);

    if (_chemical_solver_interface)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ComponentTransportLocalAssemblerInterface::setChemicalSystemID,
            _local_assemblers, pv.getActiveElementIDs());

        _chemical_solver_interface->initialize();
    }

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &ComponentTransportLocalAssemblerInterface::getIntPtDarcyVelocity));
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

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::initializeChemicalSystem,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t);
}

void ComponentTransportProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
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
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b);
}

void ComponentTransportProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ComponentTransportProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        process_id, M, K, b, Jac);
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

void ComponentTransportProcess::
    setCoupledTermForTheStaggeredSchemeToLocalAssemblers(int const process_id)
{
    DBUG("Set the coupled term for the staggered scheme to local assemblers.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ComponentTransportLocalAssemblerInterface::
            setStaggeredCoupledSolutions,
        _local_assemblers, pv.getActiveElementIDs(), _coupled_solutions);
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
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    if (process_id == 0)
    {
        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ComponentTransportLocalAssemblerInterface::setChemicalSystem,
            _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t, dt);

        BaseLib::RunTime time_phreeqc;
        time_phreeqc.start();

        _chemical_solver_interface->setAqueousSolutionsPrevFromDumpFile();

        _chemical_solver_interface->executeSpeciationCalculation(dt);

        INFO("[time] Phreeqc took {:g} s.", time_phreeqc.elapsed());

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ComponentTransportLocalAssemblerInterface::
                postSpeciationCalculation,
            _local_assemblers, pv.getActiveElementIDs(), t, dt);

        return;
    }

    auto const matrix_specification =
        ode_sys.getMatrixSpecifications(process_id);

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
    GlobalVector const& x_dot,
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
        x_dot, process_id);
}

void ComponentTransportProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
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
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, x, t, dt);

    if (!_surfaceflux)  // computing the surfaceflux is optional
    {
        return;
    }
    _surfaceflux->integrate(x, t, *this, process_id, _integration_order, _mesh,
                            pv.getActiveElementIDs());
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
