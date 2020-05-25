/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessVariable.h"

#include <algorithm>
#include <utility>
#include "BaseLib/Logging.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/TimeInterval.h"
#include "MeshGeoToolsLib/ConstructMeshesFromGeometries.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/CreateBoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/DirichletBoundaryConditionWithinTimeInterval.h"
#include "ProcessLib/SourceTerms/CreateSourceTerm.h"
#include "ProcessLib/SourceTerms/SourceTerm.h"

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
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
    :  //! \ogs_file_param{prj__process_variables__process_variable__name}
      name_(config.getConfigParameter<std::string>("name")),
      mesh_(mesh),
      //! \ogs_file_param{prj__process_variables__process_variable__components}
      n_components_(config.getConfigParameter<int>("components")),
      //! \ogs_file_param{prj__process_variables__process_variable__order}
      shapefunction_order_(config.getConfigParameter<unsigned>("order")),
      deactivated_subdomains_(createDeactivatedSubdomains(config, mesh)),
      initial_condition_(ParameterLib::findParameter<double>(
          //! \ogs_file_param{prj__process_variables__process_variable__initial_condition}
          config.getConfigParameter<std::string>("initial_condition"),
          parameters, n_components_, &mesh))
{
    DBUG("Constructing process variable {:s}", name_);

    if (shapefunction_order_ < 1 || 2 < shapefunction_order_)
    {
        OGS_FATAL("The given shape function order {:d} is not supported",
                  shapefunction_order_);
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
            auto const& mesh = findMeshInConfig(bc_config, meshes);
            auto component_id =
                //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__component}
                bc_config.getConfigParameterOptional<int>("component");

            if (!component_id && n_components_ == 1)
            {
                // default value for single component vars.
                component_id = 0;
            }

            bc_configs_.emplace_back(std::move(bc_config), mesh, component_id);
        }
    }
    else
    {
        INFO("No boundary conditions for process variable '{:s}' found.",
             name_);
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
            MeshLib::Mesh const& mesh = findMeshInConfig(st_config, meshes);
            auto component_id =
                //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__component}
                st_config.getConfigParameterOptional<int>("component");

            if (!component_id && n_components_ == 1)
            {
                // default value for single component vars.
                component_id = 0;
            }

            source_term_configs_.emplace_back(std::move(st_config), mesh,
                                              component_id);
        }
    }
    else
    {
        INFO("No source terms for process variable '{:s}' found.", name_);
    }
}

ProcessVariable::ProcessVariable(ProcessVariable&& other)
    : name_(std::move(other.name_)),
      mesh_(other.mesh_),
      n_components_(other.n_components_),
      shapefunction_order_(other.shapefunction_order_),
      deactivated_subdomains_(std::move(other.deactivated_subdomains_)),
      initial_condition_(std::move(other.initial_condition_)),
      bc_configs_(std::move(other.bc_configs_)),
      source_term_configs_(std::move(other.source_term_configs_))
{
}

std::string const& ProcessVariable::getName() const
{
    return name_;
}

MeshLib::Mesh const& ProcessVariable::getMesh() const
{
    return mesh_;
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
    bcs.reserve(bc_configs_.size());

    for (auto& config : bc_configs_)
    {
        auto bc = createBoundaryCondition(
            config, dof_table, mesh_, variable_id, integration_order,
            shapefunction_order_, parameters, process);
#ifdef USE_PETSC
        if (bc == nullptr)
        {
            continue;
        }
#endif  // USE_PETSC
        bcs.push_back(std::move(bc));
    }

    if (deactivated_subdomains_.empty())
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
    auto& parameter = ParameterLib::findParameter<double>(
        DeactivatedSubdomain::zero_parameter_name, parameters, 1);

    for (auto const& deactivated_subdomain : deactivated_subdomains_)
    {
        auto const& deactivated_subdomain_meshes =
            deactivated_subdomain->deactivated_subdomain_meshes;
        for (auto const& deactivated_subdomain_mesh :
             deactivated_subdomain_meshes)
        {
            for (int component_id = 0;
                 component_id < dof_table.getNumberOfComponents();
                 component_id++)
            {
                // Copy the time interval.
                std::unique_ptr<BaseLib::TimeInterval> time_interval =
                    std::make_unique<BaseLib::TimeInterval>(
                        *deactivated_subdomain->time_interval);

                auto bc = std::make_unique<
                    DirichletBoundaryConditionWithinTimeInterval>(
                    std::move(time_interval), parameter,
                    *(deactivated_subdomain_mesh->mesh),
                    deactivated_subdomain_mesh->inactive_nodes, dof_table,
                    variable_id, component_id);

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
    if (deactivated_subdomains_.empty())
    {
        ids_of_active_elements_.clear();
        return;
    }

    auto found_a_set =
        std::find_if(deactivated_subdomains_.begin(),
                     deactivated_subdomains_.end(),
                     [&](auto& deactivated_subdomain_) {
                         return deactivated_subdomain_->includesTimeOf(time);
                     });

    if (found_a_set == deactivated_subdomains_.end())
    {
        ids_of_active_elements_.clear();
        return;
    }

    // Already initialized.
    if (!ids_of_active_elements_.empty())
    {
        return;
    }

    auto const& deactivated_materialIDs = (*found_a_set)->materialIDs;

    auto const* const material_ids = MeshLib::materialIDs(mesh_);
    ids_of_active_elements_.clear();
    auto const number_of_elements = mesh_.getNumberOfElements();

    for (std::size_t i = 0; i < number_of_elements; i++)
    {
        if (std::binary_search(deactivated_materialIDs.begin(),
                               deactivated_materialIDs.end(),
                               (*material_ids)[i]))
        {
            continue;
        }
        ids_of_active_elements_.push_back(mesh_.getElement(i)->getID());
    }
}

std::vector<std::unique_ptr<SourceTerm>> ProcessVariable::createSourceTerms(
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const int variable_id,
    unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    std::vector<std::unique_ptr<SourceTerm>> source_terms;

    for (auto& config : source_term_configs_)
    {
        source_terms.emplace_back(createSourceTerm(
            config, dof_table, config.mesh, variable_id, integration_order,
            shapefunction_order_, parameters));
    }

    return source_terms;
}

}  // namespace ProcessLib
