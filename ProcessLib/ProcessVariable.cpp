/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessVariable.h"

#include <utility>
#include <logog/include/logog.hpp>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
ProcessVariable::ProcessVariable(
    BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh,
    GeoLib::GEOObjects const& geometries,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
    :  //! \ogs_file_param{prj__process_variables__process_variable__name}
      _name(config.getConfigParameter<std::string>("name")),
      _mesh(mesh),
      //! \ogs_file_param{prj__process_variables__process_variable__components}
      _n_components(config.getConfigParameter<int>("components")),
      //! \ogs_file_param{prj__process_variables__process_variable__order}
      _shapefunction_order(config.getConfigParameter<unsigned>("order")),
      _initial_condition(findParameter<double>(
          //! \ogs_file_param{prj__process_variables__process_variable__initial_condition}
          config.getConfigParameter<std::string>("initial_condition"),
          parameters, _n_components)),
      _bc_builder(new BoundaryConditionBuilder())
{
    DBUG("Constructing process variable %s", _name.c_str());

    if (_shapefunction_order < 1 || 2 < _shapefunction_order)
        OGS_FATAL("The given shape function order %d is not supported", _shapefunction_order);

    // Boundary conditions
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions}
    if (auto bcs_config = config.getConfigSubtreeOptional("boundary_conditions"))
    {
        for (auto bc_config :
             //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition}
             bcs_config->getConfigSubtreeList("boundary_condition"))
        {
            auto const geometrical_set_name =
                    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__geometrical_set}
                    bc_config.getConfigParameter<std::string>("geometrical_set");
            auto const geometry_name =
                    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__geometry}
                    bc_config.getConfigParameter<std::string>("geometry");

            GeoLib::GeoObject const* const geometry =
                geometries.getGeoObject(geometrical_set_name, geometry_name);

            if (! geometry)
                OGS_FATAL(
                    "No geometry with name `%s' has been found in the "
                    "geometrical set `%s'.",
                    geometry_name.c_str(), geometrical_set_name.c_str());

            DBUG(
                "Found geometry type \"%s\"",
                GeoLib::convertGeoTypeToString(geometry->getGeoType()).c_str());

            auto component_id =
                //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__component}
                bc_config.getConfigParameterOptional<int>("component");

            if (!component_id)
            {
                if (_n_components == 1)
                    // default value for single component vars.
                    component_id = 0;
                else
                    OGS_FATAL(
                        "The <component> tag could not be found for the "
                        "multi-component boundary condition of the process "
                        "variable `%s'.",
                        _name.c_str());
            }

            _bc_configs.emplace_back(std::move(bc_config), *geometry,
                                     *component_id);
        }
    } else {
        INFO("No boundary conditions found.");
    }
}

ProcessVariable::ProcessVariable(ProcessVariable&& other)
    : _name(std::move(other._name)),
      _mesh(other._mesh),
      _n_components(other._n_components),
      _shapefunction_order(other._shapefunction_order),
      _initial_condition(std::move(other._initial_condition)),
      _bc_configs(std::move(other._bc_configs)),
      _bc_builder(std::move(other._bc_builder))
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
    if (_mesh.getProperties().hasPropertyVector(_name))
    {
        auto result =
            _mesh.getProperties().template getPropertyVector<double>(_name);
        assert(result);
        assert(result->size() == _mesh.getNumberOfNodes() * _n_components);
        return *result;
    }
    else
    {
        auto result = _mesh.getProperties().template createNewPropertyVector<double>(
            _name, MeshLib::MeshItemType::Node, _n_components);
        assert(result);
        result->resize(_mesh.getNumberOfNodes() * _n_components);
        return *result;
    }
}

std::vector<std::unique_ptr<BoundaryCondition>>
ProcessVariable::createBoundaryConditions(
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const int variable_id,
    unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    std::vector<std::unique_ptr<BoundaryCondition>> bcs;

    for (auto& config : _bc_configs)
        bcs.emplace_back(_bc_builder->createBoundaryCondition(
            config, dof_table, _mesh, variable_id, integration_order, _shapefunction_order, parameters));

    return bcs;
}

}  // namespace ProcessLib
