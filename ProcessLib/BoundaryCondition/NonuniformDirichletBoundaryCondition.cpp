/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonuniformDirichletBoundaryCondition.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
std::unique_ptr<NonuniformDirichletBoundaryCondition>
createNonuniformDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, MeshLib::Mesh const& bulk_mesh)
{
    DBUG("Constructing NonuniformDirichlet BC from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NonuniformDirichlet");

    // The axial symmetry is not used in the Dirichlet BC but kept here for
    // consistency.
    //
    // Surface mesh and bulk mesh must have equal axial symmetry flags!
    if (boundary_mesh.isAxiallySymmetric() != bulk_mesh.isAxiallySymmetric())
    {
        OGS_FATAL(
            "The boundary mesh %s axially symmetric but the bulk mesh %s. Both "
            "must have an equal axial symmetry property.",
            boundary_mesh.isAxiallySymmetric() ? "is" : "is not",
            bulk_mesh.isAxiallySymmetric() ? "is" : "is not");
    }

    // TODO finally use ProcessLib::Parameter here
    auto const field_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformDirichlet__field_name}
        config.getConfigParameter<std::string>("field_name");

    auto const* const property =
        boundary_mesh.getProperties().getPropertyVector<double>(field_name);

    if (!property)
    {
        OGS_FATAL("A property with name `%s' does not exist in `%s'.",
                  field_name.c_str(), boundary_mesh.getName().c_str());
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

    return std::make_unique<NonuniformDirichletBoundaryCondition>(
        boundary_mesh, *property, dof_table, variable_id, component_id);
}

}  // ProcessLib
