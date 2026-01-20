// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ProcessVariable.h"

#include <algorithm>
#include <utility>

#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "MeshGeoToolsLib/ConstructMeshesFromGeometries.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/BoundaryCondition.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/BoundaryConditionConfig.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/CreateBoundaryCondition.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/CreateSourceTerm.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/DeactivatedSubdomainDirichlet.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/SourceTerm.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/SourceTermConfig.h"
#include "ProcessLib/CreateDeactivatedSubdomain.h"
#include "ProcessLib/DeactivatedSubdomain.h"

namespace
{
std::vector<std::reference_wrapper<const MeshLib::Mesh>> findMeshInConfig(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes)
{
    //
    // Get the mesh name from the config.
    //
    std::vector<std::string>
        mesh_names;  // Either given directly in <mesh> or <meshes> or
                     // constructed from <geometrical_set>_<geometry>.

#ifdef DOXYGEN_DOCU_ONLY
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__mesh}
    config.getConfigParameterOptional<std::string>("mesh");
#endif  // DOXYGEN_DOCU_ONLY

    auto optional_mesh_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__mesh}
        config.getConfigParameterOptional<std::string>("mesh");

    if (optional_mesh_name)
    {
        mesh_names.emplace_back(*optional_mesh_name);
        INFO("Configure mesh '{}' for boundary condition or source term.",
             mesh_names.back());
    }
    else if (
        auto const geometrical_set_name =
            //! \ogs_file_param_special{prj__process_variables__process_variable__source_terms__source_term__geometrical_set}
            //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__geometrical_set}
        config.getConfigParameterOptional<std::string>("geometrical_set"))
    {
#ifdef DOXYGEN_DOCU_ONLY
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__geometrical_set}
        config.getConfigParameterOptional<std::string>("geometrical_set");
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__geometry}
        config.getConfigParameter<std::string>("geometry");
#endif  // DOXYGEN_DOCU_ONLY

        auto const geometry_name =
            //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__geometry}
            config.getConfigParameter<std::string>("geometry");

        mesh_names.emplace_back(MeshGeoToolsLib::meshNameFromGeometry(
            *geometrical_set_name, geometry_name));
    }
    else if (
        auto const meshes_config =
            //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__meshes}
        config.getConfigSubtreeOptional("meshes"))
    {
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__meshes__mesh}
        for (auto mesh_config : meshes_config->getConfigParameterList("mesh"))
        {
            mesh_names.push_back(mesh_config.getValue<std::string>());
            INFO("Configure mesh '{:s}' for boundary condition.",
                 mesh_names.back());
        }
    }

    //
    // Find and extract mesh from the list of meshes.
    //
    std::vector<std::reference_wrapper<const MeshLib::Mesh>> bc_meshes;
    for (auto const& mesh_name : mesh_names)
    {
        bc_meshes.push_back(MeshLib::findMeshByName(meshes, mesh_name));
        DBUG("Found mesh '{:s}' with id {:d}.", mesh_name,
             bc_meshes.back().get().getID());
    }

    return bc_meshes;
}
}  // namespace

namespace ProcessLib
{
ProcessVariable::ProcessVariable(
    BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh,
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
          parameters, _n_components, &mesh)),
      _compensate_non_equilibrium_initial_residuum(
          //! \ogs_file_param{prj__process_variables__process_variable__compensate_non_equilibrium_initial_residuum}
          config.getConfigParameter<bool>(
              "compensate_non_equilibrium_initial_residuum", false))
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
            auto const bc_meshes = findMeshInConfig(bc_config, meshes);
            auto component_id =
                //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__component}
                bc_config.getConfigParameterOptional<int>("component");

            if (!component_id && _n_components == 1)
            {
                // default value for single component vars.
                component_id = 0;
            }

            _bc_configs.emplace_back(
                std::move(bc_config), bc_meshes, component_id,
                _compensate_non_equilibrium_initial_residuum);
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
            auto st_meshes = findMeshInConfig(st_config, meshes);
            auto component_id =
                //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__component}
                st_config.getConfigParameterOptional<int>("component");
            auto const type =
                //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
                st_config.peekConfigParameter<std::string>("type");

            if (!component_id)
            {
                if (_n_components == 1)
                {
                    // default value for single component vars.
                    component_id = 0;
                }
                else if (type == "Anchor")
                {
                    // dummy value
                    component_id = 0;
                }
                else if (type == "EmbeddedAnchor")
                {
                    // dummy value
                    component_id = 0;
                }
                else
                {
                    OGS_FATAL(
                        "Specifying the component id (<component>) for a "
                        "source term for a non-scalar process variable is "
                        "mandatory.");
                }
            }

            _source_term_configs.emplace_back(std::move(st_config),
                                              st_meshes.front(), *component_id);
        }
    }
    else
    {
        INFO("No source terms for process variable '{:s}' found.", _name);
    }

    if (!_deactivated_subdomains.empty())
    {
        _is_active = getOrCreateMeshProperty<unsigned char>(
            _mesh, _name + "_active", MeshLib::MeshItemType::Cell, 1);
        std::fill(std::begin(*_is_active), std::end(*_is_active), 1u);
    }
}

ProcessVariable::ProcessVariable(ProcessVariable&&) = default;

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
    Process const& process,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::vector<std::unique_ptr<BoundaryCondition>> bcs;
    bcs.reserve(_bc_configs.size());

    for (auto const& config : _bc_configs)
    {
        auto bc_vec = createBoundaryCondition(
            config, dof_table, _mesh, variable_id, integration_order,
            _shapefunction_order, parameters, process,
            all_process_variables_for_this_process, media);
        for (auto& bc : bc_vec)
        {
            // bc can be a nullptr if in parallel execution the bc isn't
            // associated to the particular subdomain
            if (bc != nullptr)
            {
                bcs.push_back(std::move(bc));
            }
        }
    }

    createBoundaryConditionsForDeactivatedSubDomains(dof_table, variable_id,
                                                     parameters, bcs);

    // Clear configs to trigger ConfigTree destruction and validation of
    // unprocessed parameters before the simulation starts.
    _bc_configs.clear();

    return bcs;
}

void ProcessVariable::createBoundaryConditionsForDeactivatedSubDomains(
    const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::vector<std::unique_ptr<BoundaryCondition>>& bcs)
{
    for (auto const& deactivated_subdomain : _deactivated_subdomains)
    {
        auto const& deactivated_subdomain_mesh =
            deactivated_subdomain.deactivated_subdomain_mesh;
        auto const* parameter = &ParameterLib::findParameter<double>(
            DeactivatedSubdomain::zero_parameter_name, parameters, 1);
        bool const set_outer_nodes_dirichlet_values =
            deactivated_subdomain.boundary_value_parameter != nullptr;
        if (set_outer_nodes_dirichlet_values)
        {
            parameter = deactivated_subdomain.boundary_value_parameter;
        }

        for (int component_id = 0;
             component_id <
             dof_table.getNumberOfVariableComponents(variable_id);
             component_id++)
        {
            auto bc = std::make_unique<DeactivatedSubdomainDirichlet>(
                *_is_active, deactivated_subdomain.time_interval, *parameter,
                set_outer_nodes_dirichlet_values, deactivated_subdomain_mesh,
                dof_table, variable_id, component_id);

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

void ProcessVariable::updateDeactivatedSubdomains(double const time)
{
    if (_deactivated_subdomains.empty())
    {
        return;
    }

    _ids_of_active_elements.clear();

    // If none of the deactivated subdomains is active at current time, then the
    // _ids_of_active_elements remain empty.
    if (std::none_of(begin(_deactivated_subdomains),
                     end(_deactivated_subdomains), [&](auto const& ds)
                     { return ds.isInTimeSupportInterval(time); }))
    {
        // Also mark all of the elements as active.
        assert(_is_active != nullptr);  // guaranteed by constructor
        std::fill(std::begin(*_is_active), std::end(*_is_active), 1u);

        return;
    }

    auto is_active_in_subdomain = [&](std::size_t const element_id,
                                      DeactivatedSubdomain const& ds) -> bool
    {
        return (!ds.isInTimeSupportInterval(time)) ||
               !ds.isDeactivated(*_mesh.getElement(element_id), time);
    };

    auto is_active_in_all_subdomains = [&](std::size_t const element_id) -> bool
    {
        return std::all_of(begin(_deactivated_subdomains),
                           end(_deactivated_subdomains),
                           [&](auto const& ds)
                           { return is_active_in_subdomain(element_id, ds); });
    };

    auto const number_of_elements = _mesh.getNumberOfElements();
    for (std::size_t element_id = 0; element_id < number_of_elements;
         element_id++)
    {
        if (is_active_in_all_subdomains(element_id))
        {
            _ids_of_active_elements.push_back(element_id);
        }
    }

    // all elements are deactivated
    std::fill(std::begin(*_is_active), std::end(*_is_active), 0u);

    for (auto const id : _ids_of_active_elements)
    {
        (*_is_active)[id] = 1u;
    }
}

std::vector<std::unique_ptr<SourceTermBase>> ProcessVariable::createSourceTerms(
    const NumLib::LocalToGlobalIndexMap& dof_table, const int variable_id,
    unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process,
    const MeshLib::Mesh& bulk_mesh)
{
    std::vector<std::unique_ptr<SourceTermBase>> source_terms;

    transform(cbegin(_source_term_configs), cend(_source_term_configs),
              back_inserter(source_terms),
              [&](auto const& config)
              {
                  return createSourceTerm(
                      config, dof_table, config.mesh, variable_id,
                      integration_order, _shapefunction_order, parameters,
                      all_process_variables_for_this_process, bulk_mesh);
              });

    // Clear configs to trigger ConfigTree destruction and validation of
    // unprocessed parameters before the simulation starts.
    _source_term_configs.clear();

    return source_terms;
}

ProcessVariable::~ProcessVariable() = default;

}  // namespace ProcessLib
