/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TESPROCESS_H_
#define PROCESS_LIB_TESPROCESS_H_

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
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
template <typename GlobalSetup>
class TESProcess final : public Process<GlobalSetup>
{
    using BP = Process<GlobalSetup>;  //!< "Base Process"

public:
    using GlobalVector = typename GlobalSetup::VectorType;
    using GlobalMatrix = typename GlobalSetup::MatrixType;

    TESProcess(
        MeshLib::Mesh& mesh,
        typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&&
            time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SecondaryVariableCollection<GlobalVector>&& secondary_variables,
        ProcessOutput<GlobalVector>&& process_output,
        BaseLib::ConfigTree const& config);

    void preTimestep(GlobalVector const& x, const double t,
                     const double delta_t) override;
    void preIteration(const unsigned iter, GlobalVector const& x) override;
    NumLib::IterationResult postIteration(GlobalVector const& x) override;

    bool isLinear() const override { return false; }
private:
    using LocalAssembler =
        TESLocalAssemblerInterface<GlobalMatrix, GlobalVector>;

    using GlobalAssembler = NumLib::VectorMatrixAssembler<
        GlobalMatrix, GlobalVector, LocalAssembler,
        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    using ExtrapolatorInterface =
        NumLib::Extrapolator<GlobalVector, TESIntPtVariables, LocalAssembler>;
    using ExtrapolatorImplementation =
        NumLib::LocalLinearLeastSquaresExtrapolator<
            GlobalVector, TESIntPtVariables, LocalAssembler>;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    GlobalVector const& computeVapourPartialPressure(
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache);

    GlobalVector const& computeRelativeHumidity(
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache);

    GlobalVector const& computeEquilibriumLoading(
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache);

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;

    AssemblyParams _assembly_params;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    std::unique_ptr<ExtrapolatorInterface> _extrapolator;

    // used for checkBounds()
    std::unique_ptr<GlobalVector> _x_previous_timestep;
};

template <typename GlobalSetup>
std::unique_ptr<TESProcess<GlobalSetup>> createTESProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&&
        time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
    BaseLib::ConfigTree const& config)
{
    config.checkConfParam("type", "TES");

    DBUG("Create TESProcess.");

    auto process_variables = findProcessVariables(
        variables, config,
        {"fluid_pressure", "temperature", "vapour_mass_fraction"});

    SecondaryVariableCollection<typename GlobalSetup::VectorType>
        secondary_variables{
            config.getConfSubtreeOptional("secondary_variables"),
            {"solid_density", "reaction_rate", "velocity_x", "velocity_y",
             "velocity_z", "loading", "reaction_damping_factor",
             "vapour_partial_pressure", "relative_humidity",
             "equilibrium_loading"}};

    ProcessOutput<typename GlobalSetup::VectorType> process_output{
        config.getConfSubtree("output"), process_variables,
        secondary_variables};

    return std::unique_ptr<TESProcess<GlobalSetup>>{new TESProcess<GlobalSetup>{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(process_variables), std::move(secondary_variables),
        std::move(process_output), config}};
}

}  // namespace TES

}  // namespace ProcessLib

#endif  // PROCESS_LIB_TESPROCESS_H_
