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

#include "ProcessLib/LIE/SmallDeformation/LocalAssembler/CreateLocalAssemblers.h"
#include "ProcessLib/LIE/SmallDeformation/LocalAssembler/SmallDeformationLocalAssemblerInterface.h"

#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{

template <int DisplacementDim>
class SmallDeformationProcess final : public Process
{
    using Base = Process;

    static_assert(DisplacementDim==2,
                  "Currently LIE::SmallDeformationProcess "
                  "supports only 2D.");

public:
    SmallDeformationProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    ProcessType getProcessType() const override
                     {return ProcessLib::ProcessType::SmallDeformationProcess;}

private:
    using LocalAssemblerInterface = SmallDeformationLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b,
                                 StaggeredCouplingTerm const&
                                 coupled_term) override
    {
        DBUG("Assemble SmallDeformationProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b,
            coupled_term);
    }

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupled_term) override
    {
        DBUG("AssembleWithJacobian SmallDeformationProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, coupled_term);
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                     double const dt) override
    {
        DBUG("PreTimestep SmallDeformationProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &SmallDeformationLocalAssemblerInterface::preTimestep,
            _local_assemblers, *_local_to_global_index_map, x, t, dt);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override;

private:
    SmallDeformationProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    std::vector<MeshLib::Element*> _vec_matrix_elements;
    std::vector<int> _vec_fracture_mat_IDs;
    std::vector<std::vector<MeshLib::Element*>> _vec_fracture_elements;
    std::vector<std::vector<MeshLib::Element*>> _vec_fracture_matrix_elements;
    std::vector<std::vector<MeshLib::Node*>> _vec_fracture_nodes;

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>> _mesh_subset_fracture_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_matrix_nodes;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
