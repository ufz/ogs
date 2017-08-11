/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonuniformNeumannBoundaryCondition.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
std::unique_ptr<NonuniformNeumannBoundaryCondition>
createNonuniformNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, MeshLib::Mesh const& bulk_mesh)
{
    DBUG("Constructing NonuniformNeumann BC from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NonuniformNeumann");

    // TODO handle paths correctly
    auto const mesh_file = config.getConfigParameter<std::string>("mesh");

    std::unique_ptr<MeshLib::Mesh> boundary_mesh(
        MeshLib::IO::readMeshFromFile(mesh_file));

    if (!boundary_mesh)
    {
        OGS_FATAL("Error reading mesh `%s'", mesh_file.c_str());
    }

    // Surface mesh and bulk mesh must have equal axial symmetry flags!
    boundary_mesh->setAxiallySymmetric(bulk_mesh.isAxiallySymmetric());

    // TODO finally use ProcessLib::Parameter here
    auto const field_name =
        config.getConfigParameter<std::string>("field_name");

    auto const* const property =
        boundary_mesh->getProperties().getPropertyVector<double>(field_name);

    if (!property)
    {
        OGS_FATAL("A property with name `%s' does not exist in `%s'.",
                  field_name.c_str(), mesh_file.c_str());
    }

    if (property->getMeshItemType() != MeshLib::MeshItemType::Node)
    {
        OGS_FATAL(
            "Only nodal fields are supported for non-uniform fields. Field "
            "`%s' is not nodal.",
            field_name.c_str());
    }

    if (property->getNumberOfComponents() != 1)
    {
        OGS_FATAL("`%s' is not a one-component field.", field_name.c_str());
    }

    std::string mapping_to_bulk_nodes_property = "OriginalSubsurfaceNodeIDs";
    auto const* const mapping_to_bulk_nodes =
        boundary_mesh->getProperties().getPropertyVector<std::size_t>(
            mapping_to_bulk_nodes_property);

    if (!(mapping_to_bulk_nodes &&
          mapping_to_bulk_nodes->getMeshItemType() ==
              MeshLib::MeshItemType::Node) &&
        mapping_to_bulk_nodes->getNumberOfComponents() == 1)
    {
        OGS_FATAL("Field `%s' is not set up properly.",
                  mapping_to_bulk_nodes_property.c_str());
    }

    return std::make_unique<NonuniformNeumannBoundaryCondition>(
        integration_order, shapefunction_order, bulk_mesh.getDimension(),
        std::move(boundary_mesh),
        NonuniformNeumannBoundaryConditionData{
            *property, bulk_mesh.getID(), *mapping_to_bulk_nodes, dof_table,
            variable_id, component_id});
}

}  // ProcessLib
