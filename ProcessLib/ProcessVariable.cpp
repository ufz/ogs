/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
      _bc_builder(std::make_unique<BoundaryConditionBuilder>()),
      _source_term_builder(std::make_unique<SourceTermBuilder>())
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

            if (!component_id && _n_components == 1)
                // default value for single component vars.
                component_id = 0;

            _bc_configs.emplace_back(std::move(bc_config), *geometry,
                                     component_id);
        }
    } else {
        INFO("No boundary conditions found.");
    }

    // Source terms
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms}
    if (auto sts_config = config.getConfigSubtreeOptional("source_terms"))
    {
        for (auto st_config :
             //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term}
             sts_config->getConfigSubtreeList("source_term"))
        {
            auto const geometrical_set_name =
                    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__geometrical_set}
                   st_config.getConfigParameter<std::string>("geometrical_set");
            auto const geometry_name =
                    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__geometry}
                    st_config.getConfigParameter<std::string>("geometry");

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
                //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__component}
                st_config.getConfigParameterOptional<int>("component");

            if (!component_id && _n_components == 1)
                // default value for single component vars.
                component_id = 0;

            _source_term_configs.emplace_back(std::move(st_config), *geometry,
                                              component_id);
        }
    } else {
        INFO("No source terms found.");
    }
}

ProcessVariable::ProcessVariable(ProcessVariable&& other)
    : _name(std::move(other._name)),
      _mesh(other._mesh),
      _n_components(other._n_components),
      _shapefunction_order(other._shapefunction_order),
      _initial_condition(std::move(other._initial_condition)),
      _bc_configs(std::move(other._bc_configs)),
      _bc_builder(std::move(other._bc_builder)),
      _source_term_configs(std::move(other._source_term_configs)),
      _source_term_builder(std::move(other._source_term_builder))
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
    return *MeshLib::getOrCreateMeshProperty<double>(
        _mesh, _name, MeshLib::MeshItemType::Node, _n_components);
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
            config, dof_table, _mesh, variable_id, integration_order,
            _shapefunction_order, parameters));

    return bcs;
}

std::vector<std::unique_ptr<NodalSourceTerm>>
ProcessVariable::createSourceTerms(
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const int variable_id,
    unsigned const integration_order)
{
    std::vector<std::unique_ptr<NodalSourceTerm>> source_terms;

    for (auto& config : _source_term_configs)
        source_terms.emplace_back(_source_term_builder->createSourceTerm(
            config, dof_table, _mesh, variable_id, integration_order,
            _shapefunction_order));

    return source_terms;
}

}  // namespace ProcessLib
