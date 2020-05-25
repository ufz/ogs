/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    static_assert(DisplacementDim == 2 || DisplacementDim == 3,
                  "Currently LIE::SmallDeformationProcess "
                  "supports only 2D or 3D.");

public:
    SmallDeformationProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          GlobalVector const& x,
                                          GlobalVector const& x_dot,
                                          int const process_id) override;

private:
    using LocalAssemblerInterface = SmallDeformationLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    const int process_id) override;

private:
    SmallDeformationProcessData<DisplacementDim> process_data_;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> local_assemblers_;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        local_to_global_index_map_single_component_;

    std::vector<MeshLib::Element*> vec_matrix_elements_;
    std::vector<int> vec_fracture_mat_IDs_;
    std::vector<std::vector<MeshLib::Element*>> vec_fracture_elements_;
    std::vector<std::vector<MeshLib::Element*>> vec_fracture_matrix_elements_;
    std::vector<std::vector<MeshLib::Node*>> vec_fracture_nodes_;
    std::vector<MeshLib::Node*> vec_junction_nodes_;
    std::vector<std::vector<MeshLib::Element*>>
        vec_junction_fracture_matrix_elements_;

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>>
        mesh_subset_fracture_nodes_;
    std::unique_ptr<MeshLib::MeshSubset const> mesh_subset_junction_nodes_;
    std::unique_ptr<MeshLib::MeshSubset const> mesh_subset_matrix_nodes_;
};

extern template class SmallDeformationProcess<2>;
extern template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
