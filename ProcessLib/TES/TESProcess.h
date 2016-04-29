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

#include "ProcessLib/Process.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"

#include "TESAssemblyParams.h"
#include "TESLocalAssembler.h"

namespace MeshLib
{
    class Element;
    class Mesh;
    template <typename PROP_VAL_TYPE> class PropertyVector;
}

namespace ProcessLib
{

namespace TES
{

template<typename GlobalSetup>
class TESProcess final
        : public Process<GlobalSetup>
{
    using BP = Process<GlobalSetup>; //!< "Base Process"

public:
    using GlobalVector = typename GlobalSetup::VectorType;
    using GlobalMatrix = typename GlobalSetup::MatrixType;

    TESProcess(MeshLib::Mesh& mesh,
               typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
               std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
               std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
               SecondaryVariableCollection<GlobalVector>&& secondary_variables,
               ProcessOutput<GlobalVector>&& process_output,
               BaseLib::ConfigTree const& config);

    void preTimestep(GlobalVector const& x, const double t, const double delta_t) override;
    void preIteration(const unsigned iter, GlobalVector const& x) override;
    NumLib::IterationResult postIteration(GlobalVector const& x) override;

    bool isLinear() const override { return false; }

private:
    using LocalAssembler  = TESLocalAssemblerInterface<GlobalMatrix, GlobalVector>;

    using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector, LocalAssembler,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    using ExtrapolatorInterface = NumLib::Extrapolator<
        GlobalVector, TESIntPtVariables, LocalAssembler>;
    using ExtrapolatorImplementation = NumLib::LocalLinearLeastSquaresExtrapolator<
        GlobalVector, TESIntPtVariables, LocalAssembler>;


    // TODO move body to cpp
    void initializeConcreteProcess(
            AssemblerLib::LocalToGlobalIndexMap const& dof_table,
            MeshLib::Mesh const& mesh, unsigned const integration_order) override
    {
        DBUG("Create global assembler.");
        _global_assembler.reset(new GlobalAssembler(dof_table));

        ProcessLib::createLocalAssemblers<GlobalSetup, TESLocalAssembler>(
                    mesh.getDimension(), mesh.getElements(),
                    dof_table, integration_order, _local_assemblers,
                    _assembly_params);

        // TODO move the two data members somewhere else.
        // for extrapolation of secondary variables
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
                    new MeshLib::MeshSubsets(BP::_mesh_subset_all_nodes.get()));
        _local_to_global_index_map_single_component.reset(
                    new AssemblerLib::LocalToGlobalIndexMap(
                        std::move(all_mesh_subsets_single_component),
                        // by location order is needed for output
                        AssemblerLib::ComponentOrder::BY_LOCATION)
                    );

        _extrapolator.reset(new ExtrapolatorImplementation(
            { 0u, 0u, nullptr, _local_to_global_index_map_single_component.get(), &mesh }));

        // secondary variables
        auto add2nd = [&](
            std::string const& var_name, unsigned const n_components,
            SecondaryVariableFunctions<GlobalVector>&& fcts)
        {
            BP::_secondary_variables.addSecondaryVariable(
                        var_name, n_components, std::move(fcts));
        };
        auto makeEx = [&](TESIntPtVariables var)
        {
            return ProcessLib::makeExtrapolator(var, *_extrapolator, _local_assemblers);
        };

        add2nd("solid_density",  1, makeEx(TESIntPtVariables::SOLID_DENSITY));
        add2nd("reaction_rate",  1, makeEx(TESIntPtVariables::REACTION_RATE));
        add2nd("velocity_x",     1, makeEx(TESIntPtVariables::VELOCITY_X));
        if (mesh.getDimension() >= 2)
            add2nd("velocity_y", 1, makeEx(TESIntPtVariables::VELOCITY_Y));
        if (mesh.getDimension() >= 3)
            add2nd("velocity_z", 1, makeEx(TESIntPtVariables::VELOCITY_Z));

        add2nd("loading",        1, makeEx(TESIntPtVariables::LOADING));
        add2nd("reaction_damping_factor",
                                 1, makeEx(TESIntPtVariables::REACTION_DAMPING_FACTOR));

        namespace PH = std::placeholders;
        using Self = TESProcess<GlobalSetup>;

        add2nd("vapour_partial_pressure", 1,
            {std::bind(&Self::computeVapourPartialPressure, this, PH::_1, PH::_2), nullptr});
        add2nd("relative_humidity",       1,
            {std::bind(&Self::computeRelativeHumidity,      this, PH::_1, PH::_2), nullptr});
        add2nd("equilibrium_loading",     1,
            {std::bind(&Self::computeEquilibriumLoading,    this, PH::_1, PH::_2), nullptr});
    }

    void assembleConcreteProcess(
            const double t, GlobalVector const& x,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override;

    GlobalVector computeVapourPartialPressure(
            GlobalVector const& x, AssemblerLib::LocalToGlobalIndexMap const& dof_table);

    GlobalVector computeRelativeHumidity(
            GlobalVector const& x, AssemblerLib::LocalToGlobalIndexMap const& dof_table);

    GlobalVector computeEquilibriumLoading(
            GlobalVector const& x, AssemblerLib::LocalToGlobalIndexMap const& dof_table);

    SecondaryVariableFunctions<GlobalVector>
    makeExtrapolator(TESIntPtVariables const var) const;


    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;

    AssemblyParams _assembly_params;

    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map_single_component;

    //! Output global matrix/rhs after first iteration.
    std::size_t _timestep = 0;
    std::size_t _total_iteration = 0;

    std::unique_ptr<ExtrapolatorInterface> _extrapolator;
};

template <typename GlobalSetup>
std::unique_ptr<TESProcess<GlobalSetup>>
createTESProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
    BaseLib::ConfigTree const& config)
{
    config.checkConfParam("type", "TES");

    DBUG("Create TESProcess.");

    auto process_variables =
        findProcessVariables(variables, config,
            { "fluid_pressure", "temperature", "vapour_mass_fraction" });

    SecondaryVariableCollection<typename GlobalSetup::VectorType>
        secondary_variables{config.getConfSubtreeOptional("secondary_variables"),
            { "solid_density", "reaction_rate",
              "velocity_x", "velocity_y", "velocity_z",
              "loading", "reaction_damping_factor",
              "vapour_partial_pressure", "relative_humidity",
              "equilibrium_loading"
            }};

    ProcessOutput<typename GlobalSetup::VectorType>
        process_output{config.getConfSubtree("output"),
                process_variables, secondary_variables};

    return std::unique_ptr<TESProcess<GlobalSetup>>{
        new TESProcess<GlobalSetup>{
            mesh, nonlinear_solver, std::move(time_discretization),
            std::move(process_variables),
            std::move(secondary_variables),
            std::move(process_output),
            config
    }};
}

} // namespace TES

} // namespace ProcessLib

#endif  // PROCESS_LIB_TESPROCESS_H_
