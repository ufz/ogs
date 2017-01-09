/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "ProcessLib/HydroMechanics/CreateLocalAssemblers.h"
#include "ProcessLib/Process.h"

#include "HydroMechanicsFEM.h"
#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
/// Linear kinematics poro-mechanical/biphasic (fluid-solid mixture) model.
///
/// The mixture momentum balance and the mixture mass balance are solved under
/// fully saturated conditions.
template <int DisplacementDim>
class HydroMechanicsProcess final : public Process
{
    using Base = Process;

public:
    HydroMechanicsProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        HydroMechanicsProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller)
        : Process(mesh, std::move(jacobian_assembler), parameters,
                  integration_order, std::move(process_variables),
                  std::move(secondary_variables),
                  std::move(named_function_caller)),
          _process_data(std::move(process_data))
    {
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
    void constructDofTable() override
    {
        // Create single component dof in every of the mesh's nodes.
        _mesh_subset_all_nodes.reset(
            new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));
        // Create single component dof in the mesh's base nodes.
        _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
        _mesh_subset_base_nodes.reset(
            new MeshLib::MeshSubset(_mesh, &_base_nodes));

        // Collect the mesh subsets in a vector.

        // For pressure, which is the first
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
        all_mesh_subsets.push_back(std::unique_ptr<MeshLib::MeshSubsets>{
            new MeshLib::MeshSubsets{_mesh_subset_base_nodes.get()}});

        // For displacement.
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables()[1].get().getNumberOfComponents(),
            [&]() {
                return std::unique_ptr<MeshLib::MeshSubsets>{
                    new MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}};
            });

        std::vector<unsigned> const vec_n_components{1, DisplacementDim};
        _local_to_global_index_map.reset(new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION));
    }

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        ProcessLib::HydroMechanics::createLocalAssemblers<DisplacementDim,
                                                          LocalAssemblerData>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            // use displacment process variable for shapefunction order
            getProcessVariables()[1].get().getShapeFunctionOrder(),
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override
    {
        DBUG("Assemble HydroMechanicsProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b);
    }

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override
    {
        DBUG("AssembleJacobian HydroMechanicsProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac);
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt) override
    {
        DBUG("PreTimestep HydroMechanicsProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, _local_assemblers,
            *_local_to_global_index_map, x, t, dt);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        DBUG("PostTimestep HydroMechanicsProcess.");

        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::postTimestep, _local_assemblers,
            *_local_to_global_index_map, x);
    }

private:
    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;
    HydroMechanicsProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib
