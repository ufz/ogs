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

#include <set>
#include <tuple>

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
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

    ~TESProcess();

    bool isLinear() const override { return false; }

    void createLocalAssemblers() override;

private:
    void assembleConcreteProcess(
            const double t, GlobalVector const& x,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override;

    template <unsigned GlobalDim>
    void createLocalAssemblers();

    void output(std::string const& file_name,
                GlobalVector const& x);

    using LocalAssembler = TESLocalAssemblerInterface<GlobalMatrix, GlobalVector>;
    std::vector<LocalAssembler*> _local_assemblers;

    AssemblerLib::ComponentOrder _global_matrix_order =
            AssemblerLib::ComponentOrder::BY_COMPONENT;

    AssemblyParams _assembly_params;

    std::vector<std::tuple<SecondaryVariables, std::string, unsigned> >
    _secondary_process_vars;

    std::set<std::string> _output_variables;

    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map_single_component;

    bool _output_residuals = false;
    //! Output global matrix/rhs after first iteration.
    bool _output_global_matrix = false;
    bool _output_iteration_results = false;

    std::size_t _timestep = 0;
    std::size_t _total_iteration = 0;

    //! Extrapolator Interface
    using ExtrapolatorIntf = NumLib::Extrapolator<GlobalVector, SecondaryVariables, LocalAssembler>;
    using ExtrapolatorImpl = NumLib::LocalLinearLeastSquaresExtrapolator<GlobalVector, SecondaryVariables, LocalAssembler>;
    std::unique_ptr<ExtrapolatorIntf> _extrapolator;
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
