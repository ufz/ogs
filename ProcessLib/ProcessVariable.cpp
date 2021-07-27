/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessVariable.h"

#include <algorithm>
#include <utility>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "MeshGeoToolsLib/ConstructMeshesFromGeometries.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/BoundaryCondition.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/CreateBoundaryCondition.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/CreateSourceTerm.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/DeactivatedSubdomainDirichlet.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/SourceTerm.h"
#include "ProcessLib/CreateDeactivatedSubdomain.h"

namespace
{
MeshLib::Mesh const& findMeshInConfig(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    //
    // Get the mesh name from the config.
    //
    std::string mesh_name;  // Either given directly in <mesh> or constructed
                            // from <geometrical_set>_<geometry>.

#ifdef DOXYGEN_DOCU_ONLY
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__mesh}
    config.getConfigParameterOptional<std::string>("mesh");
#endif  // DOXYGEN_DOCU_ONLY

    auto optional_mesh_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__mesh}
        config.getConfigParameterOptional<std::string>("mesh");
    if (optional_mesh_name)
    {
        mesh_name = *optional_mesh_name;
    }
    else
    {
#ifdef DOXYGEN_DOCU_ONLY
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__geometrical_set}
        config.getConfigParameterOptional<std::string>("geometrical_set");
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__geometry}
        config.getConfigParameter<std::string>("geometry");
#endif  // DOXYGEN_DOCU_ONLY

        // Looking for the mesh created before for the given geometry.
        auto const geometrical_set_name =
            //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__geometrical_set}
            config.getConfigParameter<std::string>("geometrical_set");
        auto const geometry_name =
            //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__geometry}
            config.getConfigParameter<std::string>("geometry");

        mesh_name = MeshGeoToolsLib::meshNameFromGeometry(geometrical_set_name,
                                                          geometry_name);
    }

    //
    // Find and extract mesh from the list of meshes.
    //
    auto const& mesh = *BaseLib::findElementOrError(
        begin(meshes), end(meshes),
        [&mesh_name](auto const& mesh) {
            assert(mesh != nullptr);
            return mesh->getName() == mesh_name;
        },
        "Required mesh with name '" + mesh_name + "' not found.");
    DBUG("Found mesh '{:s}' with id {:d}.", mesh.getName(), mesh.getID());

    return mesh;
}
}  // namespace

namespace ProcessLib
{
ProcessVariable::ProcessVariable(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh& mesh,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
    :  //! \ogs_file_param{prj__process_variables__process_variable__name}
      _name(config.getConfigParameter<std::string>("name")),
      _mesh(mesh),
      //! \ogs_file_param{prj__process_variables__process_variable__components}
      _n_components(config.getConfigParameter<int>("components")),
      //! \ogs_file_param{prj__process_variables__process_variable__order}
      _shapefunction_order(config.getConfigParameter<unsigned>("order")),
      _deactivated_subdomains(
          createDeactivatedSubdomains(config, mesh, parameters, curves)),
      _initial_condition(ParameterLib::findParameter<double>(
          //! \ogs_file_param{prj__process_variables__process_variable__initial_condition}
          config.getConfigParameter<std::string>("initial_condition"),
          parameters, _n_components, &mesh))
{
    DBUG("Constructing process variable {:s}", _name);

    if (_shapefunction_order < 1 || 2 < _shapefunction_order)
    {
        OGS_FATAL("The given shape function order {:d} is not supported",
                  _shapefunction_order);
    }

    // Boundary conditions
    if (auto bcs_config =
            //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions}
        config.getConfigSubtreeOptional("boundary_conditions"))
    {
        for (
            auto bc_config :
            //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition}
            bcs_config->getConfigSubtreeList("boundary_condition"))
        {
            auto const& bc_mesh = findMeshInConfig(bc_config, meshes);
            auto component_id =
                //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__component}
                bc_config.getConfigParameterOptional<int>("component");

            if (!component_id && _n_components == 1)
            {
                // default value for single component vars.
                component_id = 0;
            }

            _bc_configs.emplace_back(std::move(bc_config), bc_mesh,
                                     component_id);
        }
    }
    else
    {
        INFO("No boundary conditions for process variable '{:s}' found.",
             _name);
    }

    // Source terms
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms}
    if (auto sts_config = config.getConfigSubtreeOptional("source_terms"))
    {
        for (
            auto st_config :
            //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term}
            sts_config->getConfigSubtreeList("source_term"))
        {
            MeshLib::Mesh const& st_mesh = findMeshInConfig(st_config, meshes);
            auto component_id =
                //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__component}
                st_config.getConfigParameterOptional<int>("component");

            if (!component_id && _n_components == 1)
            {
                // default value for single component vars.
                component_id = 0;
            }

            _source_term_configs.emplace_back(std::move(st_config), st_mesh,
                                              component_id);
        }
    }
    else
    {
        INFO("No source terms for process variable '{:s}' found.", _name);
    }
}

ProcessVariable::ProcessVariable(ProcessVariable&& other)
    : _name(std::move(other._name)),
      _mesh(other._mesh),
      _n_components(other._n_components),
      _shapefunction_order(other._shapefunction_order),
      _deactivated_subdomains(std::move(other._deactivated_subdomains)),
      _initial_condition(std::move(other._initial_condition)),
      _bc_configs(std::move(other._bc_configs)),
      _source_term_configs(std::move(other._source_term_configs))
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

std::vector<std::unique_ptr<BoundaryCondition>>
ProcessVariable::createBoundaryConditions(
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const int variable_id,
    unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    Process const& process)
{
    std::vector<std::unique_ptr<BoundaryCondition>> bcs;
    bcs.reserve(_bc_configs.size());

    for (auto& config : _bc_configs)
    {
        auto bc = createBoundaryCondition(
            config, dof_table, _mesh, variable_id, integration_order,
            _shapefunction_order, parameters, process);
#ifdef USE_PETSC
        if (bc == nullptr)
        {
            continue;
        }
#endif  // USE_PETSC
        bcs.push_back(std::move(bc));
    }

    if (_deactivated_subdomains.empty())
    {
        return bcs;
    }

    createBoundaryConditionsForDeactivatedSubDomains(dof_table, variable_id,
                                                     parameters, bcs);
    return bcs;
}

void ProcessVariable::createBoundaryConditionsForDeactivatedSubDomains(
    const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::vector<std::unique_ptr<BoundaryCondition>>& bcs)
{
    for (auto const& deactivated_subdomain : _deactivated_subdomains)
    {
        auto const& deactivated_subdomain_meshes =
            deactivated_subdomain->deactivated_subdomain_meshes;
        for (auto const& deactivated_subdomain_mesh :
             deactivated_subdomain_meshes)
        {
            auto const* parameter = &ParameterLib::findParameter<double>(
                DeactivatedSubdomain::zero_parameter_name, parameters, 1);
            bool const set_outer_nodes_dirichlet_values =
                deactivated_subdomain->boundary_value_parameter != nullptr;
            if (set_outer_nodes_dirichlet_values)
            {
                parameter = deactivated_subdomain->boundary_value_parameter;
            }

            for (int component_id = 0;
                 component_id <
                 dof_table.getNumberOfVariableComponents(variable_id);
                 component_id++)
            {
                auto bc = std::make_unique<DeactivatedSubdomainDirichlet>(
                    &_ids_of_active_elements,
                    deactivated_subdomain->time_interval, *parameter,
                    set_outer_nodes_dirichlet_values,
                    *deactivated_subdomain_mesh, dof_table, variable_id,
                    component_id);

#ifdef USE_PETSC
                // TODO: make it work under PETSc too.
                if (bc == nullptr)
                {
                    continue;
                }
#endif  // USE_PETSC
                bcs.push_back(std::move(bc));
            }
        }
    }
}

void ProcessVariable::updateDeactivatedSubdomains(double const time)
{
    _ids_of_active_elements.clear();

    // If none of the deactivated subdomains is active at current time, then the
    // _ids_of_active_elements remain empty.
    if (std::none_of(
            begin(_deactivated_subdomains), end(_deactivated_subdomains),
            [&](auto const& ds) { return ds->isInTimeSupportInterval(time); }))
    {
        return;
    }

    auto const* const material_ids = MeshLib::materialIDs(_mesh);

    auto is_active_in_subdomain = [&](std::size_t const i,
                                      DeactivatedSubdomain const& ds) -> bool {
        if (!ds.isInTimeSupportInterval(time))
        {
            return true;
        }

        auto const& deactivated_materialIDs = ds.materialIDs;

        auto const& element_center = getCenterOfGravity(*_mesh.getElement(i));
        if (std::binary_search(deactivated_materialIDs.begin(),
                               deactivated_materialIDs.end(),
                               (*material_ids)[i]) &&
            ds.isDeactivated(element_center, time))
        {
            return false;
        }
        return true;
    };

    auto const number_of_elements = _mesh.getNumberOfElements();
    for (std::size_t i = 0; i < number_of_elements; i++)
    {
        if (std::all_of(
                begin(_deactivated_subdomains), end(_deactivated_subdomains),
                [&](auto const& ds) { return is_active_in_subdomain(i, *ds); }))
        {
            _ids_of_active_elements.push_back(_mesh.getElement(i)->getID());
        }
    }
}

std::vector<std::unique_ptr<SourceTerm>> ProcessVariable::createSourceTerms(
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const int variable_id,
    unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    std::vector<std::unique_ptr<SourceTerm>> source_terms;

    transform(cbegin(_source_term_configs), cend(_source_term_configs),
              back_inserter(source_terms), [&](auto const& config) {
                  return createSourceTerm(config, dof_table, config.mesh,
                                          variable_id, integration_order,
                                          _shapefunction_order, parameters);
              });

    return source_terms;
}

}  // namespace ProcessLib
