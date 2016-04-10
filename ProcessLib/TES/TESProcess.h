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
               std::vector<ProcessVariable> const& variables,
               std::vector<std::unique_ptr<ParameterBase>> const& parameters,
               BaseLib::ConfigTree const& config);

    void preTimestep(GlobalVector const& x, const double t, const double delta_t) override;
    void preIteration(const unsigned iter, GlobalVector const& x) override;
    NumLib::IterationResult postIteration(GlobalVector const& x) override;

    bool isLinear() const override { return false; }

    void createLocalAssemblers() override;

private:
    using LocalAssembler = TESLocalAssemblerInterface<GlobalMatrix, GlobalVector>;
    using ExtrapolatorInterface = NumLib::Extrapolator<
        GlobalVector, TESIntPtVariables, LocalAssembler>;
    using ExtrapolatorImplementation = NumLib::LocalLinearLeastSquaresExtrapolator<
        GlobalVector, TESIntPtVariables, LocalAssembler>;


    void assembleConcreteProcess(
            const double t, GlobalVector const& x,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override;

    template <unsigned GlobalDim>
    void createLocalAssemblers();

    GlobalVector computeVapourPartialPressure(
            GlobalVector const& x, AssemblerLib::LocalToGlobalIndexMap const& dof_table);

    GlobalVector computeRelativeHumidity(
            GlobalVector const& x, AssemblerLib::LocalToGlobalIndexMap const& dof_table);

    GlobalVector computeEquilibriumLoading(
            GlobalVector const& x, AssemblerLib::LocalToGlobalIndexMap const& dof_table);

    SecondaryVariableFunctions<GlobalVector>
    makeExtrapolator(TESIntPtVariables const var) const;


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
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    config.checkConfParam("type", "TES");

    DBUG("Create TESProcess.");

    return std::unique_ptr<TESProcess<GlobalSetup>>{
        new TESProcess<GlobalSetup>{
            mesh, nonlinear_solver, std::move(time_discretization),
            variables, parameters,
            config
    }};
}

} // namespace TES

} // namespace ProcessLib

#endif  // PROCESS_LIB_TESPROCESS_H_
