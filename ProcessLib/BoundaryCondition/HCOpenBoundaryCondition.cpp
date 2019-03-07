/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HCOpenBoundaryCondition.h"

#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
std::unique_ptr<HCOpenBoundaryCondition>
createHCOpenBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters, 
    unsigned const global_dim, Process const& process,
    unsigned const shapefunction_order)
{
    DBUG("Constructing open boundary for Component Transport process from config.");
//    //! \ogs_file_param{prj__processes__process__type}
//    config.checkConfigParameter("type", "ComponentTransport");
//    How can I check for process type == Component Transport? The method above
//    does not work.
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "HCOpenBoundary");


    auto const boundary_permeability_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__VariableDependentNeumann__constant_name}
        config.getConfigParameter<std::string>("boundary_permeability");
    auto const& boundary_permeability = findParameter<double>(boundary_permeability_name, parameters, 1);

    auto const bulk_element_ids =
        bc_mesh.getProperties().template getPropertyVector<std::size_t>(
            "bulk_element_ids", MeshLib::MeshItemType::Cell, 1);
    auto const bulk_face_ids =
        bc_mesh.getProperties().template getPropertyVector<std::size_t>(
            "bulk_face_ids", MeshLib::MeshItemType::Cell, 1);
    
    // In case of partitioned mesh the boundary could be empty, i.e. there is no
    // boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createBoundaryCondition().
    if (bc_mesh.getDimension() == 0 && bc_mesh.getNumberOfNodes() == 0 &&
        bc_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<HCOpenBoundaryCondition>(
        integration_order, shapefunction_order, dof_table, variable_id,
        component_id, global_dim, bc_mesh,
        HCOpenBoundaryConditionData{
            boundary_permeability, *bulk_face_ids, *bulk_element_ids, process});
}

}  // namespace ProcessLib
