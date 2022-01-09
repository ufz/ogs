/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "StokesFlowFEM.h"
#include "StokesFlowProcessData.h"

#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace StokesFlow
{
/// The Stokes equations are reduced from the Navier-Stokes equations,
/// and include momentum equations and a continuity equation describing
/// the flow in the free-flow region.
template <int GlobalDim>
class StokesFlowProcess final : public Process
{
public:
    StokesFlowProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        StokesFlowProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme);

    bool isLinear() const override { return false; }

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int /*process_id*/) const override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int /*process_id*/) const override;

    void computeSecondaryVariableConcrete(double const /*t*/,
                                          double const /*dt*/,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& /*x_dot*/,
                                          int const /*process_id*/) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t,
                                     const double dt,
                                     int const process_id) override;

private:
    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void initializeBoundaryConditions() override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;

    StokesFlowProcessData _process_data;

    std::vector<std::unique_ptr<StokesFlowLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    /// Local to global index mapping for base nodes, which is used for linear
    /// interpolation for pressure in the staggered scheme.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_with_base_nodes;
};

extern template class StokesFlowProcess<2>;
}  // namespace StokesFlow
}  // namespace ProcessLib
