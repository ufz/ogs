/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_H_
#define PROCESS_LIB_PROCESS_H_

#include <memory>
#include <string>

#include <logog/include/logog.hpp>

#include "AssemblerLib/ComputeSparsityPattern.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#include "BaseLib/ConfigTree.h"
#include "FileIO/VtkIO/VtuInterface.h"
#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSubsets.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#include "MathLib/LinAlg/PETSc/PETScMatrixOption.h"
#endif

#include "DirichletBc.h"
#include "NeumannBc.h"
#include "NeumannBcAssembler.h"
#include "Parameter.h"
#include "ProcessVariable.h"
#include "UniformDirichletBoundaryCondition.h"

#include "NumLib/ODESolver/ODESystem.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
template <typename GlobalSetup>
class Process
		: public NumLib::ODESystem<typename GlobalSetup::MatrixType,
		                           typename GlobalSetup::VectorType,
		                           // TODO: later on use a simpler ODE system
		                           NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
		                           NumLib::NonlinearSolverTag::Picard>
{
public:
	using GlobalVector = typename GlobalSetup::VectorType;
	using GlobalMatrix = typename GlobalSetup::MatrixType;

	Process(MeshLib::Mesh& mesh) : _mesh(mesh) {}
	virtual ~Process()
	{
		for (auto p : _all_mesh_subsets)
			delete p;
		delete _mesh_subset_all_nodes;
	}

	/// Process specific initialization called by initialize().
	virtual void createLocalAssemblers() = 0;
	virtual bool assemble(const double delta_t) = 0;

	virtual std::string getLinearSolverName() const = 0;

	/// Postprocessing after solve().
	/// The file_name is indicating the name of possible output file.
	void postTimestep(std::string const& file_name,
	                  const unsigned /*timestep*/,
	                  GlobalVector const& x)
	{
		post(x);
		output(file_name, x);
	}

	void initialize()
	{
		DBUG("Initialize process.");

		DBUG("Construct dof mappings.");
		initializeMeshSubsets();

		_local_to_global_index_map.reset(
		    new AssemblerLib::LocalToGlobalIndexMap(
		        _all_mesh_subsets, AssemblerLib::ComponentOrder::BY_COMPONENT));

#ifndef USE_PETSC
		DBUG("Compute sparsity pattern");
		computeSparsityPattern();
#endif

		// create global vectors and linear solver
		createLinearSolver(getLinearSolverName());

		DBUG("Create global assembler.");
		_global_assembler.reset(
		    new GlobalAssembler(*_local_to_global_index_map));

		createLocalAssemblers();

		DBUG("Initialize boundary conditions.");
		for (ProcessVariable& pv : _process_variables)
		{
			createDirichletBcs(pv, 0);  // 0 is the component id
			createNeumannBcs(pv, 0);    // 0 is the component id
		}

		for (auto& bc : _neumann_bcs)
			bc->initialize(_global_setup, *_A, *_rhs, _mesh.getDimension());
	}

	void setInitialConditions(GlobalVector& x)
	{
		DBUG("Set initial conditions.");
		for (ProcessVariable& pv : _process_variables)
		{
			setInitialConditions(pv, 0, x);  // 0 is the component id
		}
	}

	bool solve_TODO_DELETE(const double delta_t)
	{
		_A->setZero();
		MathLib::setMatrixSparsity(*_A, _sparsity_pattern);

		bool const result = assemble(delta_t);

		// Call global assembler for each Neumann boundary local assembler.
		for (auto const& bc : _neumann_bcs)
			bc->integrate(_global_setup);

		for (auto const& bc : _dirichlet_bcs)
			MathLib::applyKnownSolution(*_A, *_rhs, *_x, bc.global_ids,
			                            bc.values);

		_linear_solver->solve(*_rhs, *_x);
		return result;
	}

protected:
	virtual void post(GlobalVector const& x)
	{
		(void) x; // by default do nothing
	}

	/// Set linear solver options; called by the derived process which is
	/// parsing the configuration.
	void setLinearSolverOptions(BaseLib::ConfigTree&& config)
	{
		_linear_solver_options.reset(
		    new BaseLib::ConfigTree(std::move(config)));
	}

private:
	/// Creates mesh subsets, i.e. components, for given mesh.
	void initializeMeshSubsets()
	{
		// Create single component dof in every of the mesh's nodes.
		_mesh_subset_all_nodes =
		    new MeshLib::MeshSubset(_mesh, &_mesh.getNodes());

		// Collect the mesh subsets in a vector.
		_all_mesh_subsets.push_back(
		    new MeshLib::MeshSubsets(_mesh_subset_all_nodes));
	}

	/// Sets the initial condition values in the solution vector x for a given
	/// process variable and component.
	void setInitialConditions(ProcessVariable const& variable,
	                          int const component_id,
	                          GlobalVector& x)
	{
		std::size_t const n = _mesh.getNNodes();
		for (std::size_t i = 0; i < n; ++i)
		{
			MeshLib::Location const l(_mesh.getID(),
			                          MeshLib::MeshItemType::Node, i);
			auto global_index = std::abs(
			    _local_to_global_index_map->getGlobalIndex(l, component_id));
#ifdef USE_PETSC
			// The global indices of the ghost entries of the global
			// matrix or the global vectors need to be set as negative values
			// for equation assembly, however the global indices start from zero.
			// Therefore, any ghost entry with zero index is assigned an negative
			// value of the vector size or the matrix dimension.
			// To assign the initial value for the ghost entries,
			// the negative indices of the ghost entries are restored to zero.
			// checked hereby.
			if ( global_index == x.size() )
			    global_index = 0;
#endif
			x.set(global_index,
			      variable.getInitialConditionValue(*_mesh.getNode(i)));
		}
	}

	void createDirichletBcs(ProcessVariable& variable, int const component_id)
	{
		MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
		    MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
		        variable.getMesh());

		variable.initializeDirichletBCs(std::back_inserter(_dirichlet_bcs),
		                                mesh_node_searcher,
		                                *_local_to_global_index_map,
		                                component_id);
	}

	void createNeumannBcs(ProcessVariable& variable, int const component_id)
	{
		// Find mesh nodes.
		MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
		    MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
		        variable.getMesh());
		MeshGeoToolsLib::BoundaryElementsSearcher mesh_element_searcher(
		    variable.getMesh(), mesh_node_searcher);

		// Create a neumann BC for the process variable storing them in the
		// _neumann_bcs vector.
		variable.createNeumannBcs(std::back_inserter(_neumann_bcs),
		                          mesh_element_searcher,
		                          _global_setup,
		                          _integration_order,
		                          *_local_to_global_index_map,
		                          component_id,
		                          *_mesh_subset_all_nodes);
	}

	/// Creates global matrix, rhs and solution vectors, and the linear solver.
	void createLinearSolver(std::string const& solver_name)
	{
		DBUG("Allocate global matrix, vectors, and linear solver.");
#ifdef USE_PETSC
		MathLib::PETScMatrixOption mat_opt;
		const MeshLib::NodePartitionedMesh& pmesh =
		    static_cast<const MeshLib::NodePartitionedMesh&>(_mesh);
		mat_opt.d_nz = pmesh.getMaximumNConnectedNodesToNode();
		mat_opt.o_nz = mat_opt.d_nz;
		mat_opt.is_global_size = false;
		const std::size_t num_unknowns =
		    _local_to_global_index_map->dofSizeLocal();
		_A.reset(_global_setup.createMatrix(num_unknowns, mat_opt));
		// In the following two lines, false is assigned to
		// the argument of is_global_size, which indicates num_unknowns
		// is local.
		_x.reset( _global_setup.createVector(num_unknowns,
		          _local_to_global_index_map->getGhostIndices(), false) );
		_rhs.reset( _global_setup.createVector(num_unknowns,
		            _local_to_global_index_map->getGhostIndices(), false) );
#else
		const std::size_t num_unknowns = _local_to_global_index_map->dofSize();
		_A.reset(_global_setup.createMatrix(num_unknowns));
		_x.reset(_global_setup.createVector(num_unknowns));
		_rhs.reset(_global_setup.createVector(num_unknowns));
#endif
		_linear_solver.reset(new typename GlobalSetup::LinearSolver(
		    *_A, solver_name, _linear_solver_options.get()));
		checkAndInvalidate(_linear_solver_options);
	}

	/// Computes and stores global matrix' sparsity pattern from given
	/// DOF-table.
	void computeSparsityPattern()
	{
		_sparsity_pattern = std::move(AssemblerLib::computeSparsityPattern(
		    *_local_to_global_index_map, _mesh));
	}

	void output(std::string const& file_name, GlobalVector const& x)
	{
		DBUG("Process output.");

		std::string const property_name = "Result";

		// Get or create a property vector for results.
		boost::optional<MeshLib::PropertyVector<double>&> result;
		if (_mesh.getProperties().hasPropertyVector(property_name))
		{
			result = _mesh.getProperties().template getPropertyVector<double>(
			    property_name);
		}
		else
		{
			result =
			    _mesh.getProperties().template createNewPropertyVector<double>(
			        property_name, MeshLib::MeshItemType::Node);
#ifdef USE_PETSC
			result->resize(x.getLocalSize() + x.getGhostSize());
#else
			result->resize(x.size());
#endif
		}

		assert(result);

		// Copy result
		x.copyValues(*result);

		// Write output file
		DBUG("Writing output to \'%s\'.", file_name.c_str());
		FileIO::VtuInterface vtu_interface(&_mesh, vtkXMLWriter::Binary, true);
		vtu_interface.writeToFile(file_name);
	}

protected:
	unsigned const _integration_order = 2;

	MeshLib::Mesh& _mesh;
	MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
	std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

	GlobalSetup _global_setup;

	using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<
	        GlobalMatrix,
	        GlobalVector,
	        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

	std::unique_ptr<GlobalAssembler> _global_assembler;

	std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap>
	    _local_to_global_index_map;

	std::unique_ptr<BaseLib::ConfigTree> _linear_solver_options;
	std::unique_ptr<typename GlobalSetup::LinearSolver> _linear_solver;

	std::unique_ptr<GlobalMatrix> _A;
	std::unique_ptr<GlobalVector> _rhs;
	std::unique_ptr<GlobalVector> _x;

	AssemblerLib::SparsityPattern _sparsity_pattern;

	std::vector<DirichletBc<GlobalIndexType>> _dirichlet_bcs;
	std::vector<std::unique_ptr<NeumannBc<GlobalSetup>>> _neumann_bcs;

	/// Variables used by this process.
	std::vector<std::reference_wrapper<ProcessVariable>> _process_variables;
};

/// Find a process variable for a name given in the process configuration under
/// the tag.
/// In the process config a process variable is referenced by a name. For
/// example it will be looking for a variable named "H" in the list of process
/// variables when the tag is "hydraulic_head":
/// \code
///     <process>
///         ...
///         <hydraulic_head>H</hydraulic_head>
///     </process>
/// \endcode
/// and return a reference to that variable.
ProcessVariable& findProcessVariable(
    BaseLib::ConfigTree const& process_config, std::string const& tag,
    std::vector<ProcessVariable> const& variables);

/// Find a parameter of specific type for a name given in the process
/// configuration under the tag.
/// In the process config a parameter is referenced by a name. For example it
/// will be looking for a parameter named "K" in the list of parameters
/// when the tag is "hydraulic_conductivity":
/// \code
///     <process>
///         ...
///         <hydraulic_conductivity>K</hydraulic_conductivity>
///     </process>
/// \endcode
/// and return a reference to that parameter. Additionally it checks for the
/// type of the found parameter.
template <typename... ParameterArgs>
Parameter<ParameterArgs...>& findParameter(
    BaseLib::ConfigTree const& process_config, std::string const& tag,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
	// Find parameter name in process config.
	auto const name = process_config.getConfParam<std::string>(tag);

	// Find corresponding parameter by name.
	auto const parameter_it =
	    std::find_if(parameters.cbegin(), parameters.cend(),
	                 [&name](std::unique_ptr<ParameterBase> const& p)
	                 {
		                 return p->name == name;
	                 });

	if (parameter_it == parameters.end())
	{
		ERR(
		    "Could not find parameter '%s' in the provided parameters list for "
		    "config tag <%s>.",
		    name.c_str(), tag.c_str());
		std::abort();
	}
	DBUG("Found parameter \'%s\'.", (*parameter_it)->name.c_str());

	// Check the type correctness of the found parameter.
	auto* const parameter =
	    dynamic_cast<Parameter<ParameterArgs...>*>(parameter_it->get());
	if (!parameter)
	{
		ERR("The read parameter is of incompatible type.");
		std::abort();
	}
	return *parameter;
}

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_H_
