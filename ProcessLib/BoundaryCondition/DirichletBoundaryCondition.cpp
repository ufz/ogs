/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DirichletBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
DirichletBoundaryCondition::DirichletBoundaryCondition(
    Parameter<double> const& parameter, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id)
    : _parameter(parameter),
      _bc_mesh(bc_mesh),
      _variable_id(variable_id),
      _component_id(component_id)
{
    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, "
            "%d), maximum values: (%d, %d).",
            variable_id, component_id, dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    if (!_bc_mesh.getProperties().existsPropertyVector<std::size_t>(
            "bulk_node_ids"))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the boundary "
            "mesh '%s'.",
            _bc_mesh.getName().c_str());
    }

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    DBUG(
        "Found %d nodes for Dirichlet BCs for the variable %d and "
        "component "
        "%d",
        bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset)));
}

void DirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    return getEssentialBCValuesLocal(t, x, bc_values);
}

void DirichletBoundaryCondition::getEssentialBCValuesLocal(
    const double t, GlobalVector const& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + _bc_mesh.getNumberOfNodes());
    bc_values.values.reserve(bc_values.values.size() +
                             _bc_mesh.getNumberOfNodes());
    for (auto const* const node : _bc_mesh.getNodes())
    {
        auto const id = node->getID();
        pos.setNodeID(node->getID());
        // TODO: that might be slow, but only done once
        auto const global_index = _dof_table_boundary->getGlobalIndex(
            {_bc_mesh.getID(), MeshLib::MeshItemType::Node, id}, _variable_id,
            _component_id);
        if (global_index == NumLib::MeshComponentMap::nop)
            continue;
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply the
        // Dirichlet BC, the negative index is not accepted like other matrix or
        // vector PETSc routines. Therefore, the following if-condition is
        // applied.
        if (global_index >= 0)
        {
            bc_values.ids.emplace_back(global_index);
            bc_values.values.emplace_back(_parameter(t, pos).front());
        }
    }
}

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    DBUG("Constructing DirichletBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "Dirichlet");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);

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

    return std::make_unique<DirichletBoundaryCondition>(
        param, bc_mesh, dof_table_bulk, variable_id, component_id);
}

}  // namespace ProcessLib
