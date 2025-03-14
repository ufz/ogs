/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HydroMechanicsProcessData.h"
#include "LocalAssembler/HydroMechanicsLocalAssemblerInterface.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
class HydroMechanicsLocalAssemblerInterface;

template <int GlobalDim>
class HydroMechanicsProcess final : public Process
{
    static_assert(GlobalDim == 2 || GlobalDim == 3,
                  "Currently LIE::HydroMechanicsProcess "
                  "supports only 2D or 3D.");

public:
    HydroMechanicsProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        HydroMechanicsProcessData<GlobalDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     double const t, double const dt,
                                     int const process_id) override;

private:
    using LocalAssemblerInterface = HydroMechanicsLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac) override;
    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    int const process_id) override;

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
