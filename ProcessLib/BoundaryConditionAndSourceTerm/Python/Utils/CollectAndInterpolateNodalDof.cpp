/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CollectAndInterpolateNodalDof.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace
{
/// Collects the degrees of freedom of the passed element from the passed global
/// vector into a vector.
///
/// \note \c num_nodes can be set "arbitrarily", e.g., to only collect d.o.f. on
/// the base nodes of the passed mesh element.
void collectDofsToMatrixSingleComponentForSomeNodes(
    MeshLib::Element const& element, std::size_t const mesh_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    int const variable, int const component, unsigned const num_nodes,
    Eigen::Ref<Eigen::VectorXd> all_nodal_dof_for_this_component)
{
    bool dof_not_found = false;

    for (unsigned element_node_id = 0; element_node_id < num_nodes;
         ++element_node_id)
    {
        auto const& node = *element.getNode(element_node_id);
        auto const node_id = node.getID();
        MeshLib::Location const loc{mesh_id, MeshLib::MeshItemType::Node,
                                    node_id};
        auto const dof_idx = dof_table.getGlobalIndex(loc, variable, component);

        if (dof_idx == NumLib::MeshComponentMap::nop)
        {
            // We just skip this d.o.f. Actually we will also skip
            // all other nodes of this mesh element for this (var,
            // comp), because we assume that all linear nodes have a
            // lower node id than any higher order node.
            dof_not_found = true;
        }
        else
        {
            if (dof_not_found)
            {
                // We expect that for all mesh elements all linear
                // nodes have a lower node id than any higher order
                // node. I.e., there are no "holes" in the
                // primary_variables_mat. We rely on there being no
                // "holes" later on when interpolating the nodal
                // d.o.f. to the integration points.
                OGS_FATAL(
                    "This d.o.f. has been found in the d.o.f. "
                    "table, but before some d.o.f. has not been "
                    "found. Something has gone terribly wrong. "
                    "Some assumption in the implementation is "
                    "wrong.");
            }
            all_nodal_dof_for_this_component[element_node_id] = x[dof_idx];
        }
    }
}

}  // namespace

namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
{

Eigen::MatrixXd collectDofsToMatrix(
    MeshLib::Element const& element,
    std::size_t const mesh_id,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    GlobalVector const& x)
{
    auto const num_var = dof_table.getNumberOfVariables();
    auto const num_nodes = element.getNumberOfNodes();
    auto const num_comp_total = dof_table.getNumberOfGlobalComponents();

    Eigen::MatrixXd primary_variables_mat(num_nodes, num_comp_total);

    for (int var = 0; var < num_var; ++var)
    {
        auto const num_comp = dof_table.getNumberOfVariableComponents(var);

        for (int comp = 0; comp < num_comp; ++comp)
        {
            auto const global_component =
                dof_table.getGlobalComponent(var, comp);
            auto all_nodal_dof_for_this_component =
                primary_variables_mat.col(global_component);

            collectDofsToMatrixSingleComponentForSomeNodes(
                element, mesh_id, dof_table, x, var, comp, num_nodes,
                all_nodal_dof_for_this_component);
        }
    }

    return primary_variables_mat;
}

Eigen::VectorXd collectDofsToMatrixOnBaseNodesSingleComponent(
    MeshLib::Element const& element, std::size_t const mesh_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    int const variable, int const component)
{
    auto const num_nodes = element.getNumberOfBaseNodes();
    Eigen::VectorXd primary_variables_vec(num_nodes);

    collectDofsToMatrixSingleComponentForSomeNodes(
        element, mesh_id, dof_table, x, variable, component, num_nodes,
        primary_variables_vec);

    return primary_variables_vec;
}

}  // namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
