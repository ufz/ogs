/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HTProcess.h"

#include <cassert>

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "ProcessLib/CalculateSurfaceFlux/CalculateSurfaceFlux.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

#include "HTMaterialProperties.h"

#include "MonolithicHTFEM.h"
#include "StaggeredHTFEM.h"

namespace ProcessLib
{
namespace HT
{
HTProcess::HTProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    std::unique_ptr<HTMaterialProperties>&& material_properties,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme,
    std::unique_ptr<ProcessLib::Balance>&& balance)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _material_properties(std::move(material_properties)),
      _balance(std::move(balance))
{
}

void HTProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    // For the staggered scheme, both processes are assumed to use the same
    // element order. Therefore the order of shape function can be fetched from
    // any set of the sets of process variables of the coupled processes. Here,
    // we take the one from the first process by setting process_id = 0.
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    if (_use_monolithic_scheme)
    {
        ProcessLib::createLocalAssemblers<MonolithicHTFEM>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            pv.getShapeFunctionOrder(), _local_assemblers,
            mesh.isAxiallySymmetric(), integration_order,
            *_material_properties);
    }
    else
    {
        const int heat_transport_process_id = 0;
        const int hydraulic_process_id = 1;

        ProcessLib::createLocalAssemblers<StaggeredHTFEM>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            pv.getShapeFunctionOrder(), _local_assemblers,
            mesh.isAxiallySymmetric(), integration_order, *_material_properties,
            heat_transport_process_id, hydraulic_process_id);
    }

    _secondary_variables.addSecondaryVariable(
        "darcy_velocity",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         _local_assemblers,
                         &HTLocalAssemblerInterface::getIntPtDarcyVelocity));
}

void HTProcess::assembleConcreteProcess(const double t,
                                        GlobalVector const& x,
                                        GlobalMatrix& M,
                                        GlobalMatrix& K,
                                        GlobalVector& b)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (_use_monolithic_scheme)
    {
        DBUG("Assemble HTProcess.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        if (_coupled_solutions->process_id == 0)
        {
            DBUG(
                "Assemble the equations of heat transport process within "
                "HTProcess.");
        }
        else
        {
            DBUG(
                "Assemble the equations of single phase fully saturated "
                "fluid flow process within HTProcess.");
        }
        setCoupledSolutionsOfPreviousTimeStep();
        dof_tables.emplace_back(*_local_to_global_index_map);
        dof_tables.emplace_back(*_local_to_global_index_map);
    }

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_tables, t, x, M, K, b, _coupled_solutions);
}

void HTProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HTProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    if (!_use_monolithic_scheme)
    {
        setCoupledSolutionsOfPreviousTimeStep();
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
    }
    else
    {
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
        dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
    }

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);
}

void HTProcess::preTimestepConcreteProcess(GlobalVector const& x,
                                           const double /*t*/,
                                           const double /*delta_t*/,
                                           const int process_id)
{
    assert(process_id < 2);

    if (_use_monolithic_scheme)
    {
        return;
    }

    if (!_xs_previous_timestep[process_id])
    {
        _xs_previous_timestep[process_id] =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);
    }
    else
    {
        auto& x0 = *_xs_previous_timestep[process_id];
        MathLib::LinAlg::copy(x, x0);
    }

    auto& x0 = *_xs_previous_timestep[process_id];
    MathLib::LinAlg::setLocalAccessibleVector(x0);
}

void HTProcess::setCoupledTermForTheStaggeredSchemeToLocalAssemblers()
{
    DBUG("Set the coupled term for the staggered scheme to local assembers.");

    GlobalExecutor::executeMemberOnDereferenced(
        &HTLocalAssemblerInterface::setStaggeredCoupledSolutions,
        _local_assemblers, _coupled_solutions);
}

std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
HTProcess::getDOFTableForExtrapolatorData() const
{
    if (!_use_monolithic_scheme)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        const bool manage_storage = false;
        return std::make_tuple(_local_to_global_index_map.get(),
                               manage_storage);
    }

    // Otherwise construct a new DOF table.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};

    const bool manage_storage = true;
    return std::make_tuple(new NumLib::LocalToGlobalIndexMap(
                               std::move(all_mesh_subsets_single_component),
                               // by location order is needed for output
                               NumLib::ComponentOrder::BY_LOCATION),
                           manage_storage);
}

void HTProcess::setCoupledSolutionsOfPreviousTimeStep()
{
    const auto number_of_coupled_solutions =
        _coupled_solutions->coupled_xs.size();
    _coupled_solutions->coupled_xs_t0.clear();
    _coupled_solutions->coupled_xs_t0.reserve(number_of_coupled_solutions);
    const int process_id = _coupled_solutions->process_id;
    for (std::size_t i = 0; i < number_of_coupled_solutions; i++)
    {
        const auto& x_t0 = _xs_previous_timestep[process_id];
        if (x_t0 == nullptr)
        {
            OGS_FATAL(
                "Memory is not allocated for the global vector "
                "of the solution of the previous time step for the ."
                "staggered scheme.\n It can be done by overriding "
                "Process::preTimestepConcreteProcess"
                "(ref. HTProcess::preTimestepConcreteProcess) ");
        }

        MathLib::LinAlg::setLocalAccessibleVector(*x_t0);
        _coupled_solutions->coupled_xs_t0.emplace_back(x_t0.get());
    }
}

Eigen::Vector3d HTProcess::getFlux(std::size_t element_id,
                                   MathLib::Point3d const& p,
                                   double const t,
                                   GlobalVector const& x) const
{
    // fetch local_x from primary variable
    std::vector<GlobalIndexType> indices_cache;
    auto const r_c_indices = NumLib::getRowColumnIndices(
        element_id, *_local_to_global_index_map, indices_cache);
    std::vector<double> local_x(x.get(r_c_indices.rows));

    return _local_assemblers[element_id]->getFlux(p, t, local_x);
}

// this is almost a copy of the implemention in the GroundwaterFlow
void HTProcess::postTimestepConcreteProcess(GlobalVector const& x,
                                            const double t,
                                            const double /*delta_t*/,
                                            int const process_id)
{
    // For the monolithic scheme, process_id is always zero.
    if (_use_monolithic_scheme && process_id != 0)
    {
        OGS_FATAL(
            "The condition of process_id = 0 must be satisfied for "
            "monolithic HTProcess, which is a single process.");
    }
    if (!_use_monolithic_scheme && process_id != 1)
    {
        DBUG("This is the thermal part of the staggered HTProcess.");
        return;
    }
    if (!_balance)  // computing the balance is optional
    {
        return;
    }
    auto* const balance_pv = MeshLib::getOrCreateMeshProperty<double>(
        _balance->surface_mesh, _balance->property_vector_name,
        MeshLib::MeshItemType::Cell, 1);
    // initialise the PropertyVector pv with zero values
    std::fill(balance_pv->begin(), balance_pv->end(), 0.0);
    auto balance = ProcessLib::CalculateSurfaceFlux(
        _balance->surface_mesh,
        getProcessVariables(process_id)[0].get().getNumberOfComponents(),
        _integration_order);

    balance.integrate(
        x, *balance_pv, t, _mesh,
        [this](std::size_t const element_id, MathLib::Point3d const& pnt,
               double const t, GlobalVector const& x) {
            return getFlux(element_id, pnt, t, x);
        });
    // post: surface_mesh has scalar element property
    _balance->save(t);
}

}  // namespace HT
}  // namespace ProcessLib
