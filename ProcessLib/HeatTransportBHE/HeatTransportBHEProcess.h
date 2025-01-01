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

#include "HeatTransportBHEProcessData.h"
#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHEProcessAssemblerInterface.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
struct BHEMeshData;

class HeatTransportBHEProcess final : public Process
{
public:
    HeatTransportBHEProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        HeatTransportBHEProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        BHEMeshData&& bhe_mesh_data);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override
    {
        return _process_data._algebraic_BC_Setting._is_linear;
    }

    bool requiresNormalization() const override
    {
        // In the current setup, when using algebraic bc,
        // then normalization is always required
        return _process_data._algebraic_BC_Setting._use_algebraic_bc;
    }
    //! @}

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_prev,
                                          int const process_id) override;

private:
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

    void createBHEBoundaryConditionTopBottom(
        std::vector<std::vector<MeshLib::Node*>> const& all_bhe_nodes);
    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    const double t, const double dt,
                                    int const process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     const double t, const double dt,
                                     int const process_id) override;

    void algebraicBcConcreteProcess(const double t, double const dt,
                                    std::vector<GlobalVector*> const& x,
                                    std::vector<GlobalVector*> const& xdot,
                                    int const process_id, GlobalMatrix& M,
                                    GlobalMatrix& K, GlobalVector& b);

    NumLib::IterationResult postIterationConcreteProcess(
        GlobalVector const& x) override;

    HeatTransportBHEProcessData _process_data;

    std::vector<std::unique_ptr<HeatTransportBHELocalAssemblerInterface>>
        _local_assemblers;

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>>
        _mesh_subset_BHE_nodes;

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>>
        _mesh_subset_BHE_soil_nodes;
    // a vector of tuple structure containing the indices of BHE top nodes,
    // used only for algebraic boundary conditions
    // first object is the index of BHE
    // second and third object is the global indices of a pair of unknowns,
    // pointing to the inflow and outflow temperature
    std::vector<std::tuple<std::size_t, GlobalIndexType, GlobalIndexType>>
        _vec_top_BHE_node_indices;
    // a vector of tuple structure containing the indices of BHE bottom nodes,
    // used only for algebraic boundary conditions
    // same structure as the top node vector
    std::vector<std::tuple<std::size_t, GlobalIndexType, GlobalIndexType>>
        _vec_bottom_BHE_node_indices;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_soil_nodes;

    const BHEMeshData _bheMeshData;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
