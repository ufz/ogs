/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessVariable.h"

#include "logog/include/logog.hpp"

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

#include "UniformDirichletBoundaryCondition.h"
#include "InitialCondition.h"

namespace ProcessLib
{
ProcessVariable::ProcessVariable(BaseLib::ConfigTree const& config,
                                 MeshLib::Mesh const& mesh,
                                 GeoLib::GEOObjects const& geometries)
    : _name(config.get<std::string>("name")), _mesh(mesh)
{
	DBUG("Constructing process variable %s", this->_name.c_str());

	// Initial condition
	{
		auto const& ic_config = config.find("initial_condition");
		if (ic_config == config.not_found())
			INFO("No initial condition found.");

		std::string const type =
		    config.get<std::string>("initial_condition.type");
		if (type == "Uniform")
		{
			_initial_condition =
			    createUniformInitialCondition(ic_config->second);
		}
		if (type == "MeshProperty")
		{
			_initial_condition =
			    createMeshPropertyInitialCondition(ic_config->second, _mesh);
		}
		else
		{
			ERR("Unknown type of the initial condition.");
		}
	}

	// Boundary conditions
	{
		auto const& bcs_config = config.find("boundary_conditions");
		if (bcs_config == config.not_found())
			INFO("No boundary conditions found.");

		for (auto const& bc_iterator : bcs_config->second)
		{
			BaseLib::ConfigTree const& bc_config = bc_iterator.second;

			// Find corresponding GeoObject
			std::string const geometrical_set_name =
			    bc_config.get<std::string>("geometrical_set");
			std::string const geometry_name =
			    bc_config.get<std::string>("geometry");

			GeoLib::GeoObject const* const geometry =
			    geometries.getGeoObject(geometrical_set_name, geometry_name);
			DBUG(
			    "Found geometry type \"%s\"",
			    GeoLib::convertGeoTypeToString(geometry->getGeoType()).c_str());

			// Construct type dependent boundary condition
			std::string const type = bc_config.get<std::string>("type");

			if (type == "UniformDirichlet")
			{
				_dirichlet_bcs.emplace_back(
				    new UniformDirichletBoundaryCondition(geometry, bc_config));
			}
			else if (type == "UniformNeumann")
			{
				_neumann_bc_configs.emplace_back(
				    new NeumannBcConfig(geometry, bc_config));
			}
			else
			{
				ERR("Unknown type \'%s\' of the boundary condition.",
				    type.c_str());
			}
		}
	}
}

ProcessVariable::ProcessVariable(ProcessVariable&& other)
    : _name(std::move(other._name)),
      _mesh(other._mesh),
      _initial_condition(std::move(other._initial_condition)),
      _dirichlet_bcs(std::move(other._dirichlet_bcs)),
      _neumann_bc_configs(std::move(other._neumann_bc_configs))
{
}

std::string const& ProcessVariable::getName() const
{
	return _name;
}

MeshLib::Mesh const& ProcessVariable::getMesh() const
{
	return _mesh;
}

void ProcessVariable::initializeDirichletBCs(
    MeshGeoToolsLib::MeshNodeSearcher& searcher,
    std::vector<std::size_t>& global_ids, std::vector<double>& values)
{
	for (auto& bc : _dirichlet_bcs)
		bc->initialize(searcher, global_ids, values);
}

}  // namespace ProcessLib
