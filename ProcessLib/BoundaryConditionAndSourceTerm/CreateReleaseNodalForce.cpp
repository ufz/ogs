// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateReleaseNodalForce.h"

#include <limits>
#include <numeric>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>

#include "BaseLib/ConfigTree.h"
#include "BoundaryCondition.h"
#include "BoundaryConditionConfig.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "ParameterLib/Utils.h"
#include "ReleaseNodalForce.h"

namespace ProcessLib
{
std::string parseReleaseNodalForce(BoundaryConditionConfig const& bc_config)
{
    DBUG("Parse ReleaseNodalForce boundary condition.");

    if (!bc_config.compensate_non_equilibrium_initial_residuum)
    {
        OGS_FATAL(
            "The compensate_non_equilibrium_initial_residuum must be set for "
            "the ReleaseNodalForce boundary condition.");
    }

    auto const& config = bc_config.config;
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "ReleaseNodalForce");

    auto const decay_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__ReleaseNodalForce__time_decay_parameter}
        config.getConfigParameter<std::string>("time_decay_parameter");
    DBUG("decay parameter {:s}", decay_parameter_name);

    return decay_parameter_name;
}

std::unique_ptr<BoundaryCondition> createReleaseNodalForce(
    unsigned const global_dim, int const variable_id,
    std::string const& decay_parameter_name,
    BoundaryConditionConfig const& bc_config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Create ReleaseNodalForce boundary condition.");

    if (!bc_config.compensate_non_equilibrium_initial_residuum)
    {
        OGS_FATAL(
            "The compensate_non_equilibrium_initial_residuum must be set for "
            "the ReleaseNodalForce boundary condition.");
    }

    auto const& time_decay_parameter = ParameterLib::findParameter<double>(
        decay_parameter_name, parameters, 1, &bc_mesh);

    // In case of partitioned mesh the boundary could be empty, i.e. there
    // is no boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the
    // config parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of
    // the subtree and move the code up in createBoundaryCondition().
    if (bc_mesh.getDimension() == 0 && bc_mesh.getNumberOfNodes() == 0 &&
        bc_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    if (bc_mesh.getDimension() >= global_dim)
    {
        OGS_FATAL(
            "The dimension ({:d}) of the given boundary mesh '{:s}' is not "
            "lower than the bulk dimension ({:d}).",
            bc_mesh.getDimension(), bc_mesh.getName(), global_dim);
    }

    // Create component ids vector for the current variable.
    auto const& number_of_components =
        dof_table_bulk.getNumberOfVariableComponents(variable_id);
    std::vector<int> component_ids =
        ranges::views::iota(0, number_of_components) | ranges::to_vector;

    // BC mesh subset creation
    std::vector<MeshLib::Node*> const bc_nodes = bc_mesh.getNodes();
    DBUG("Found {:d} nodes for Natural BCs for the variable {:d}",
         bc_nodes.size(), variable_id);

    MeshLib::MeshSubset bc_mesh_subset(bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component ids.
    auto dof_table_boundary = dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, std::move(component_ids), std::move(bc_mesh_subset));

    ParameterLib::SpatialPosition pos;
    pos.setNodeID(0);
    pos.setCoordinates(MathLib::Point3d{});
    if (time_decay_parameter(0, pos).front() - 1.0 >
        std::numeric_limits<double>::epsilon())
    {
        OGS_FATAL(
            "The time_decay_parameter '{:s}' must be equal to 1.0 at t=0.0. "
            "and decrease over time to 0. ",
            decay_parameter_name);
    }

    return std::make_unique<ReleaseNodalForce>(
        variable_id, bc_mesh, dof_table_boundary, time_decay_parameter);
}

}  // namespace ProcessLib
