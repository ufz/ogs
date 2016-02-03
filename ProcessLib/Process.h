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
#include "ProcessVariable.h"
#include "UniformDirichletBoundaryCondition.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
template <typename GlobalSetup>
class Process
{
public:
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
	void postTimestep(std::string const& file_name, const unsigned /*timestep*/)
	{
		post();
		output(file_name);
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
		    new GlobalAssembler(*_A, *_rhs, *_local_to_global_index_map));

		createLocalAssemblers();

		for (auto const& pv : _process_variables)
		{
			DBUG("Set initial conditions.");
			setInitialConditions(*pv, 0);  // 0 is the component id

			DBUG("Initialize boundary conditions.");
			createDirichletBcs(*pv, 0);  // 0 is the component id

			createNeumannBcs(*pv, 0);  // 0 is the component id
		}

		for (auto& bc : _neumann_bcs)
			bc->initialize(_global_setup, *_A, *_rhs, _mesh.getDimension());
	}

	bool solve(const double delta_t)
	{
		_A->setZero();
		MathLib::setMatrixSparsity(*_A, _sparsity_pattern);

		bool const result = assemble(delta_t);

		// Call global assembler for each Neumann boundary local assembler.
		for (auto const& bc : _neumann_bcs)
			bc->integrate(_global_setup);

		for (auto const& bc : _dirichlet_bcs)
			MathLib::applyKnownSolution(*_A, *_rhs, *_x,
			                            bc.global_ids, bc.values);

		_linear_solver->solve(*_rhs, *_x);
		return result;
	}

protected:
	virtual void post() { };

	/// Set linear solver options; called by the derived process which is
	/// parsing the configuration.
	void setLinearSolverOptions(BaseLib::ConfigTree&& config)
	{
		_linear_solver_options.reset(new BaseLib::ConfigTree(
			std::move(config)));
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
	                          int const component_id)
	{
		std::size_t const n = _mesh.getNNodes();
		for (std::size_t i = 0; i < n; ++i)
		{
			MeshLib::Location const l(_mesh.getID(),
			                          MeshLib::MeshItemType::Node, i);
			auto const global_index = std::abs(
			    _local_to_global_index_map->getGlobalIndex(l, component_id));
			_x->set(global_index,
			        variable.getInitialConditionValue(*_mesh.getNode(i)));
		}
	}

	void createDirichletBcs(ProcessVariable& variable,
	                        int const component_id)
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
		const std::size_t num_unknowns =
		    _local_to_global_index_map->dofSizeGlobal();
		_A.reset(_global_setup.createMatrix(num_unknowns, mat_opt));
#else
		const std::size_t num_unknowns = _local_to_global_index_map->dofSize();
		_A.reset(_global_setup.createMatrix(num_unknowns));
#endif
		_x.reset(_global_setup.createVector(num_unknowns));
		_rhs.reset(_global_setup.createVector(num_unknowns));
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

	void output(std::string const& file_name)
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
			result->resize(_x->size());
		}
		assert(result && result->size() == _x->size());

#ifdef USE_PETSC
		std::unique_ptr<double[]> u(new double[_x->size()]);
		_x->getGlobalVector(u.get());  // get the global solution

		std::size_t const n = _mesh.getNNodes();
		for (std::size_t i = 0; i < n; ++i)
		{
			MeshLib::Location const l(_mesh.getID(),
			                          MeshLib::MeshItemType::Node, i);
			auto const global_index = std::abs(  // 0 is the component id.
			    _local_to_global_index_map->getGlobalIndex(l, 0));
			(*result)[i] = u[global_index];
		}
#else
		// Copy result
		for (std::size_t i = 0; i < _x->size(); ++i)
			(*result)[i] = (*_x)[i];
#endif

		// Write output file
		FileIO::VtuInterface vtu_interface(&_mesh, vtkXMLWriter::Binary, true);
		vtu_interface.writeToFile(file_name);
	}

protected:
	unsigned const _integration_order = 2;

	MeshLib::Mesh& _mesh;
	MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
	std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

	GlobalSetup _global_setup;

	using GlobalAssembler =
	    AssemblerLib::VectorMatrixAssembler<typename GlobalSetup::MatrixType,
	                                        typename GlobalSetup::VectorType>;

	std::unique_ptr<GlobalAssembler> _global_assembler;

	std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap>
	    _local_to_global_index_map;

	std::unique_ptr<BaseLib::ConfigTree> _linear_solver_options;
	std::unique_ptr<typename GlobalSetup::LinearSolver> _linear_solver;

	std::unique_ptr<typename GlobalSetup::MatrixType> _A;
	std::unique_ptr<typename GlobalSetup::VectorType> _rhs;
	std::unique_ptr<typename GlobalSetup::VectorType> _x;

	AssemblerLib::SparsityPattern _sparsity_pattern;

	std::vector<DirichletBc<GlobalIndexType>> _dirichlet_bcs;
	std::vector<std::unique_ptr<NeumannBc<GlobalSetup>>> _neumann_bcs;

    /// Variables used by this process.
	std::vector<ProcessVariable*> _process_variables;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_H_
