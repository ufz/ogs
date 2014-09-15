/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_PROCESS_VARIABLE_H_
#define OGS_PROCESS_VARIABLE_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

#include "OGS/BoundaryCondition.h"
#include "OGS/InitialCondition.h"

namespace OGS
{

/// A named process variable. Its properties includes the mesh, and the initial
/// and boundary conditions.
class ProcessVariable
{
    using ConfigTree = boost::property_tree::ptree;
public:
    ProcessVariable(ConfigTree const& config, MeshLib::Mesh const& mesh,
            GeoLib::GEOObjects const& geometries)
        : _name(config.get<std::string>("name")),
          _mesh(mesh)
    {
        DBUG("Constructing process variable %s", this->_name.c_str());

        // Initial condition
        {
            auto const& ic_config = config.find("initial_condition");
            if (ic_config == config.not_found())
                ERR("No initial condition found.");

            std::string const type =
                config.get<std::string>("initial_condition.type");
            if (type == "Uniform")
            {
                _initial_condition =
                    new OGS::UniformInitialCondition(ic_config->second);
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
                ERR("No boundary conditions found.");

            for (auto const& bc_iterator : bcs_config->second)
            {
                ConfigTree const& bc_config = bc_iterator.second;

                // Find corresponding GeoObject
                std::string const geometry_name =
                    bc_config.get<std::string>("geometry");
                std::string const patch_name =
                    bc_config.get<std::string>("patch");

                // TODO Currently only Polylines are supported for the boundary
                // conditions.
                GeoLib::GeoObject const* const geometry = geometries.getGeoObject(
                        geometry_name, GeoLib::GEOTYPE::POLYLINE, patch_name);


                // Construct type dependent boundary condition
                std::string const type = bc_config.get<std::string>("type");

                if (type == "Dirichlet")
                {
                    _boundary_conditions.emplace_back(
                        new OGS::DirichletBoundaryCondition(
                            geometry, bc_config));
                }
                else
                {
                    ERR("Unknown type \'%s\' of the boundary condition.",
                            type.c_str());
                }
            }

        }
    }

    ~ProcessVariable()
    {
        delete _initial_condition;

        for(auto p : _boundary_conditions)
            delete p;
    }

    std::string const& getName() const
    {
        return _name;
    }

private:
    std::string const _name;
    MeshLib::Mesh const& _mesh;
    OGS::InitialCondition* _initial_condition;
    std::vector<OGS::BoundaryCondition*> _boundary_conditions;
};

}   // namespace OGS

#endif  // OGS_PROCESS_VARIABLE_H_
