/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        SecondaryVariableCollection&& secondary_variables);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override { return false; }

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_dot,
                                          int const process_id) override;

private:
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
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void createBHEBoundaryConditionTopBottom(
        std::vector<std::vector<MeshLib::Node*>> const& all_bhe_nodes);
#ifdef OGS_USE_PYTHON
    NumLib::IterationResult postIterationConcreteProcess(
        GlobalVector const& x) override;
#endif

    HeatTransportBHEProcessData _process_data;

    std::vector<std::unique_ptr<HeatTransportBHELocalAssemblerInterface>>
        _local_assemblers;

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>>
        _mesh_subset_BHE_nodes;

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>>
        _mesh_subset_BHE_soil_nodes;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_soil_nodes;

    const BHEMeshData _bheMeshData;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
