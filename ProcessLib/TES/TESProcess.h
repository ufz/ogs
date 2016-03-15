/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TESPROCESS_H_
#define PROCESS_LIB_TESPROCESS_H_

#include <memory>
#include <vector>
#include <set>
#include <array>
#include <tuple>

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/ComputeSparsityPattern.h"

#include "FileIO/VtkIO/VtuInterface.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/Scaling.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"
#include "MathLib/Nonlinear/Picard.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#include "ProcessLib/ProcessVariable.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Parameter.h"

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/Extrapolation/GlobalLinearLeastSquaresExtrapolator.h"

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


/// Global ids in the global matrix/vector where the dirichlet bc is
/// imposed and their corresponding values.
struct DirichletBC
{
    std::vector<GlobalIndexType> global_ids;
    std::vector<double> values;
};


template<typename GlobalSetup>
class TESProcess final
        : public Process<GlobalSetup>,
          public TESProcessInterface
{
    using BP = Process<GlobalSetup>;

    unsigned const _integration_order = 2;

public:
    using GlobalVector = typename GlobalSetup::VectorType;
    using GlobalMatrix = typename GlobalSetup::MatrixType;

    TESProcess(MeshLib::Mesh& mesh,
               typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
               std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
               std::vector<ProcessVariable> const& variables,
               std::vector<std::unique_ptr<ParameterBase>> const& parameters,
               BaseLib::ConfigTree const& config);

    void post(std::string const& file_name);
    void postTimestep(std::string const& file_name, const unsigned timestep);

    ~TESProcess();

    bool isLinear() const override { return false; }

    void createLocalAssemblers() override;

private:
    void assembleConcreteProcess(
            const double t, GlobalVector const& x,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override;

    template <unsigned GlobalDim>
    void createLocalAssemblers();

    using LocalAssembler = TESLocalAssemblerInterface<GlobalMatrix, GlobalVector>;
    using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<GlobalMatrix, GlobalVector,
    NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    GlobalSetup _global_setup;
    std::vector<LocalAssembler*> _local_assemblers;
    std::unique_ptr<GlobalAssembler> _global_assembler;

    AssemblerLib::ComponentOrder _global_matrix_order =
            AssemblerLib::ComponentOrder::BY_COMPONENT;
    bool _output_residuals = false;


    std::unique_ptr<GlobalVector> _x;           // current iteration


    // secondary variables
    std::vector<std::tuple<SecondaryVariables, std::string, unsigned> >
    _secondary_process_vars;

    // output variables
    std::set<std::string> _output_variables;

    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map_single_component;

    bool _output_global_matrix = false; ///< output global matrix/rhs at first iteration
    bool _output_iteration_results = false;

    std::size_t _timestep = 0;
    std::size_t _total_iteration = 0;


#if 0
    NumLib::GlobalLinearLeastSquaresExtrapolator<
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType, SecondaryVariables,
            LocalAssembler>
            extrapolator(*_local_to_global_index_map_single_component);
#else
    using ExtrapolatorIntf = NumLib::Extrapolator<GlobalVector, SecondaryVariables, LocalAssembler>;
    using ExtrapolatorImpl = NumLib::LocalLinearLeastSquaresExtrapolator<GlobalVector, SecondaryVariables, LocalAssembler>;
    std::unique_ptr<ExtrapolatorIntf> _extrapolator;
#endif


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
