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

#include "ProcessVariable.h"
#include "Process.h"

#include "Parameter.h"

#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/Extrapolation/GlobalLinearLeastSquaresExtrapolator.h"

#include "TESProcess-notpl.h"
#include "TESFEM.h"


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
class TESProcess
        : public Process<GlobalSetup>,
          public TESProcessInterface
{
    using BP = Process<GlobalSetup>;

    unsigned const _integration_order = 2;

public:
    using GlobalVector = typename GlobalSetup::VectorType;
    using GlobalMatrix = typename GlobalSetup::MatrixType;

    TESProcess(MeshLib::Mesh& mesh,
               std::vector<ProcessVariable> const& variables,
               std::vector<std::unique_ptr<ParameterBase>> const& parameters,
               BaseLib::ConfigTreeNew const& config);

    void init() override;

    bool assemble(/*const double current_time,*/ const double delta_t) override;

    void post(std::string const& file_name) override;
    void postTimestep(std::string const& file_name, const unsigned timestep) override;

    std::string getLinearSolverName() const override { return "tes_"; }
    void initializeMeshSubsets(MeshLib::Mesh const& /*mesh*/) override {}

    ~TESProcess();

private:
    template <unsigned GlobalDim>
    void createLocalAssemblers();

    void setInitialConditions(ProcessVariable const& variable,
                              std::size_t const component_id);

    void singlePicardIteration(GlobalVector& x_prev_iter, GlobalVector& x_curr);

    using LocalAssembler = TES::LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>;
    using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<GlobalMatrix, GlobalVector>;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    GlobalSetup _global_setup;
    std::vector<LocalAssembler*> _local_assemblers;
    std::unique_ptr<GlobalAssembler> _global_assembler;

    AssemblerLib::SparsityPattern _sparsity_pattern;

    std::unique_ptr<MathLib::Nonlinear::Picard> _picard;

    AssemblerLib::ComponentOrder _global_matrix_order =
            AssemblerLib::ComponentOrder::BY_COMPONENT;
    bool _output_residuals = false;


    // primary variables
    std::array<ProcessVariable*, NODAL_DOF> _process_vars; // ) = { nullptr, nullptr, nullptr };
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    std::unique_ptr<GlobalMatrix> _A;
    std::unique_ptr<GlobalVector> _rhs;
    std::unique_ptr<GlobalVector> _x;           // current iteration
    std::unique_ptr<GlobalVector> _x_prev_ts;   // previous timestep

    std::unique_ptr<BaseLib::ConfigTreeNew> _linear_solver_options;
    std::unique_ptr<typename GlobalSetup::LinearSolver> _linear_solver;

    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    DirichletBC _dirichlet_bc;
    std::vector<NeumannBc<GlobalSetup>*> _neumann_bcs;


    // secondary variables
    std::vector<std::tuple<SecondaryVariables, std::string, unsigned> >
    _secondary_process_vars;

    // output variables
    std::set<std::string> _output_variables;

    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets_single_component;
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

} // namespace TES

} // namespace ProcessLib

#include "TESProcess-impl.h"

#endif  // PROCESS_LIB_TESPROCESS_H_
