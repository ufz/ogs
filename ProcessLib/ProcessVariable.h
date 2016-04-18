/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_VARIABLE_H_
#define PROCESS_LIB_PROCESS_VARIABLE_H_


#include "InitialCondition.h"
#include "UniformDirichletBoundaryCondition.h"
#include "NeumannBc.h"
#include "DirichletBc.h"

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;
class BoundaryElementsSearcher;
}

namespace MeshLib
{
class Mesh;
class Node;
}

namespace GeoLib
{
class GEOObjects;
}

namespace ProcessLib
{
class NeumannBcConfig;
class InitialCondition;
class UniformDirichletBoundaryCondition;
}

namespace ProcessLib
{
/// A named process variable. Its properties includes the mesh, and the initial
/// and boundary conditions.
class ProcessVariable
{
public:
	ProcessVariable(BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh,
	                GeoLib::GEOObjects const& geometries);

	ProcessVariable(ProcessVariable&&);

	std::string const& getName() const;

	/// Returns a mesh on which the process variable is defined.
	MeshLib::Mesh const& getMesh() const;

	/// Returns the number of components of the process variable.
	int getNumberOfComponents() const { return _n_components; }

	template <typename OutputIterator>
	void initializeDirichletBCs(
	    OutputIterator output_bcs,
	    MeshGeoToolsLib::MeshNodeSearcher& searcher,
	    const AssemblerLib::LocalToGlobalIndexMap& dof_table,
	    const unsigned nodal_dof_idx)
	{
		for (auto& bc_config : _dirichlet_bc_configs)
		{
			DirichletBc<GlobalIndexType> bc;
			bc_config->initialize(searcher, dof_table, nodal_dof_idx, bc);
			output_bcs++ = bc;
		}
	}

	template <typename GlobalSetup, typename OutputIterator, typename... Args>
	void createNeumannBcs(OutputIterator bcs,
	                      MeshGeoToolsLib::BoundaryElementsSearcher& searcher,
	                      Args&&... args)
	{
		for (auto& config : _neumann_bc_configs)
		{
			config->initialize(searcher);
			bcs++ = std::unique_ptr<NeumannBc<GlobalSetup>>{
			    new NeumannBc<GlobalSetup>(*config,
			                               std::forward<Args>(args)...)};
		}
	}

	double getInitialConditionValue(MeshLib::Node const& n,
	                                int const component_id) const
	{
		return _initial_condition->getValue(n, component_id);
	}

	// Get or create a property vector for results.
	// The returned mesh property size is number of mesh nodes times number of
	// components.
	MeshLib::PropertyVector<double>& getOrCreateMeshProperty();

private:
	std::string const _name;
	MeshLib::Mesh& _mesh;
	const int _n_components;
	std::unique_ptr<InitialCondition> _initial_condition;
	std::vector<std::unique_ptr<UniformDirichletBoundaryCondition>>
	    _dirichlet_bc_configs;
	std::vector<std::unique_ptr<NeumannBcConfig>> _neumann_bc_configs;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_VARIABLE_H_
