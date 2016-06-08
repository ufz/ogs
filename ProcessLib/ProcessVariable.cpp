/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessVariable.h"

#include "logog/include/logog.hpp"

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
ProcessVariable::ProcessVariable(BaseLib::ConfigTree const& config,
                                 MeshLib::Mesh& mesh,
                                 GeoLib::GEOObjects const& geometries)
    :
    //! \ogs_file_param{prj__process_variables__process_variable__name}
    _name(config.getConfParam<std::string>("name")),
    _mesh(mesh),
    //! \ogs_file_param{prj__process_variables__process_variable__components}
    _n_components(config.getConfParam<int>("components"))
{
    DBUG("Constructing process variable %s", this->_name.c_str());

    // Initial condition
    //! \ogs_file_param{prj__process_variables__process_variable__initial_condition}
    if (auto ic_config = config.getConfSubtreeOptional("initial_condition"))
    {
        //! \ogs_file_param{initial_condition__type}
        auto const type = ic_config->peekConfParam<std::string>("type");
        if (type == "Uniform")
        {
            _initial_condition =
                createUniformInitialCondition(*ic_config, _n_components);
        }
        else if (type == "MeshProperty")
        {
            _initial_condition =
                createMeshPropertyInitialCondition(*ic_config, _mesh, _n_components);
        }
        else
        {
            ERR("Unknown type of the initial condition.");
        }
    }
    else
    {
        INFO("No initial condition found.");
    }

    // Boundary conditions
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions}
    if (auto bcs_config = config.getConfSubtreeOptional("boundary_conditions"))
    {
        for (auto bc_config :
             //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition}
             bcs_config->getConfSubtreeList("boundary_condition"))
        {
            auto const geometrical_set_name =
                    //! \ogs_file_param{boundary_condition__geometrical_set}
                    bc_config.getConfParam<std::string>("geometrical_set");
            auto const geometry_name =
                    //! \ogs_file_param{boundary_condition__geometry}
                    bc_config.getConfParam<std::string>("geometry");

            GeoLib::GeoObject const* const geometry =
                geometries.getGeoObject(geometrical_set_name, geometry_name);
            DBUG(
                "Found geometry type \"%s\"",
                GeoLib::convertGeoTypeToString(geometry->getGeoType()).c_str());

            // Construct type dependent boundary condition
            //! \ogs_file_param{boundary_condition__type}
            auto const type = bc_config.peekConfParam<std::string>("type");

            if (type == "UniformDirichlet")
            {
                _dirichlet_bc_configs.emplace_back(
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
    } else {
        INFO("No boundary conditions found.");
    }

    // Source Terms
    config.ignoreConfParam("source_terms");
}

ProcessVariable::ProcessVariable(ProcessVariable&& other)
    : _name(std::move(other._name)),
      _mesh(other._mesh),
      _n_components(other._n_components),
      _initial_condition(std::move(other._initial_condition)),
      _dirichlet_bc_configs(std::move(other._dirichlet_bc_configs)),
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

MeshLib::PropertyVector<double>& ProcessVariable::getOrCreateMeshProperty()
{
    boost::optional<MeshLib::PropertyVector<double>&> result;
    if (_mesh.getProperties().hasPropertyVector(_name))
    {
        result =
            _mesh.getProperties().template getPropertyVector<double>(_name);
        assert(result);
        assert(result->size() == _mesh.getNNodes() * _n_components);
    }
    else
    {
        result = _mesh.getProperties().template createNewPropertyVector<double>(
            _name, MeshLib::MeshItemType::Node);
        assert(result);
        result->resize(_mesh.getNNodes() * _n_components);
    }
    return *result;
}

}  // namespace ProcessLib
