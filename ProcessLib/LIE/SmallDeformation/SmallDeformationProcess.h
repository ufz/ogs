/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/LIE/SmallDeformation/LocalAssembler/SmallDeformationLocalAssemblerInterface.h"
#include "ProcessLib/Process.h"
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

    static_assert(DisplacementDim == 2 || DisplacementDim == 3,
                  "Currently LIE::SmallDeformationProcess "
                  "supports only 2D or 3D.");

public:
    SmallDeformationProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

private:
    using LocalAssemblerInterface = SmallDeformationLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt,
                                    const int /*process_id*/) override;

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

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>>
        _mesh_subset_fracture_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_matrix_nodes;
};

extern template class SmallDeformationProcess<2>;
extern template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
