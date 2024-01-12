/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HeatConductionFEM.h"
#include "HeatConductionProcessData.h"
#include "ProcessLib/Assembly/AssembledMatrixCache.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace HeatConduction
{
class HeatConductionProcess final : public Process
{
public:
    HeatConductionProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        HeatConductionProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const is_linear,
        bool const ls_compute_only_upon_timestep_change);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return _asm_mat_cache.isLinear(); }
    //! @}

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_prev,
                                          int const process_id) override;

    bool shouldLinearSolverComputeOnlyUponTimestepChange() const override
    {
        // TODO move this feature to some common location for all processes.
        return _ls_compute_only_upon_timestep_change;
    }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const /*dt*/,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const /*dt*/,
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    void preOutputConcreteProcess(const double t, double const dt,
                                  std::vector<GlobalVector*> const& x,
                                  std::vector<GlobalVector*> const& x_prev,
                                  int const process_id) override;

    HeatConductionProcessData _process_data;

    std::vector<std::unique_ptr<HeatConductionLocalAssemblerInterface>>
        _local_assemblers;

    MeshLib::PropertyVector<double>* _heat_flux = nullptr;

    AssembledMatrixCache _asm_mat_cache;

    bool const _ls_compute_only_upon_timestep_change;
};

}  // namespace HeatConduction
}  // namespace ProcessLib
