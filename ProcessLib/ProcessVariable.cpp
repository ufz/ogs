/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

#include "BoundaryCondition.h"
#include "InitialCondition.h"

#include "ProcessVariable.h"

namespace ProcessLib
{

ProcessVariable::ProcessVariable(
    ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    GeoLib::GEOObjects const& geometries)
    : _name(config.get<std::string>("name")),
      _mesh(mesh)
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
                new UniformInitialCondition(ic_config->second);
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
            ConfigTree const& bc_config = bc_iterator.second;

            // Find corresponding GeoObject
            std::string const geometrical_set_name =
                bc_config.get<std::string>("geometrical_set");
            std::string const geometry_name =
                bc_config.get<std::string>("geometry");

            GeoLib::GeoObject const* const geometry = geometries.getGeoObject(
                    geometrical_set_name, geometry_name);
            DBUG("Found geometry type \"%s\"",
                GeoLib::convertGeoTypeToString(geometry->getGeoType()).c_str());

            // Construct type dependent boundary condition
            std::string const type = bc_config.get<std::string>("type");

            if (type == "UniformDirichlet")
            {
                _boundary_conditions.emplace_back(
                    new UniformDirichletBoundaryCondition(
                        geometry, bc_config));
            }
            else
            {
                ERR("Unknown type \'%s\' of the boundary condition.",
                        type.c_str());
            }

            _boundary_conditions.back()->applyToMesh(mesh);
        }

    }
}

ProcessVariable::~ProcessVariable()
{
    delete _initial_condition;

    for(auto p : _boundary_conditions)
        delete p;
}

std::string const& ProcessVariable::getName() const
{
    return _name;
}

}   // namespace ProcessLib
