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

#include "ProcessLib/Process.h"

#include "ProcessLib/LIE/HydroMechanics/LocalAssembler/CreateLocalAssemblers.h"
#include "ProcessLib/LIE/HydroMechanics/LocalAssembler/HydroMechanicsLocalAssemblerInterface.h"

#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <unsigned GlobalDim>
class HydroMechanicsProcess final : public Process
{
    using Base = Process;

    static_assert(GlobalDim == 2 || GlobalDim == 3,
                  "Currently LIE::HydroMechanicsProcess "
                  "supports only 2D or 3D.");

public:
    HydroMechanicsProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        HydroMechanicsProcessData<GlobalDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    void computeSecondaryVariableConcrete(double const t,
                                          GlobalVector const& x) override;

private:
    using LocalAssemblerInterface = HydroMechanicsLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override
    {
        DBUG("Assemble HydroMechanicsProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b,
            _coupled_solutions);
    }

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override
    {
        DBUG("AssembleWithJacobian HydroMechanicsProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, _coupled_solutions);
    }

    void preTimestepConcreteProcess(
        GlobalVector const& x, double const t, double const dt,
        const int /*process_id*/) override
    {
        DBUG("PreTimestep HydroMechanicsProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &HydroMechanicsLocalAssemblerInterface::preTimestep,
            _local_assemblers, *_local_to_global_index_map, x, t, dt);
    }

private:
    HydroMechanicsProcessData<GlobalDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::vector<MeshLib::Element*> _vec_matrix_elements;
    std::vector<MeshLib::Element*> _vec_fracture_elements;
    std::vector<MeshLib::Element*> _vec_fracture_matrix_elements;
    std::vector<MeshLib::Node*> _vec_fracture_nodes;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_fracture_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_matrix_nodes;

    std::vector<MeshLib::Node*> _mesh_nodes_p;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_nodes_p;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
