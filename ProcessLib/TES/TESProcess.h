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

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"

#include "TESAssemblyParams.h"
#include "TESLocalAssembler.h"

namespace MeshLib
{
class Element;
class Mesh;
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
namespace TES
{
class TESProcess final : public Process
{
public:
    TESProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SecondaryVariableCollection&& secondary_variables,
        BaseLib::ConfigTree const& config);

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    const double t, const double delta_t,
                                    const int process_id) override;
    void preIterationConcreteProcess(const unsigned iter,
                                     GlobalVector const& x) override;
    NumLib::IterationResult postIterationConcreteProcess(
        GlobalVector const& x) override;

    bool isLinear() const override { return false; }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void initializeSecondaryVariables();

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    GlobalVector const& computeVapourPartialPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache);

    GlobalVector const& computeRelativeHumidity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache);

    GlobalVector const& computeEquilibriumLoading(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache);

    std::vector<std::unique_ptr<TESLocalAssemblerInterface>> _local_assemblers;

    AssemblyParams _assembly_params;

    // used for checkBounds()
    std::unique_ptr<GlobalVector> _x_previous_timestep;
};

}  // namespace TES

}  // namespace ProcessLib
