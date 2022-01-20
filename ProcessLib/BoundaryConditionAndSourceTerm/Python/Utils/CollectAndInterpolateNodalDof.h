/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
{
/// Collects the degrees of freedom of the passed element from the passed global
/// vector into a matrix.
///
/// The dimensions of the returned matrix are \#nodes x \#total_components,
/// i.e., each row of the matrix will contain all degrees of freedom at a
/// specific node of the passed element.
///
/// The order of nodes is determined by the passed element. The order of
/// components is determined by the passed dof_table.
///
/// \note OGS currently has implemented Lagrange finite elements only, i.e., all
/// degrees of freedom are nodal degrees of freedom.
///
/// \note In the case of Taylor-Hood elements, some components are defined only
/// on base nodes. In that case, the returned matrix might contain some
/// uninitialized entries. It is the responsibility of the caller to handle the
/// returned matrix correctly.
Eigen::MatrixXd collectDofsToMatrix(
    MeshLib::Element const& element,
    std::size_t const mesh_id,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    GlobalVector const& x);

/// Does the same as collectDofsToMatrix(), just for a single component and only
/// on the base nodes of the passed mesh element.
Eigen::VectorXd collectDofsToMatrixOnBaseNodesSingleComponent(
    MeshLib::Element const& element, std::size_t const mesh_id,
    NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
    int const variable, int const component);

/// Interpolates the passed matrix of degrees of freedom to the \c
/// interpolated_primary_variables via the passed shape matrices in \c
/// ns_and_weight.
///
/// The interpolation function order for each variable is determined by the
/// passed ProcessVariable vector.
///
/// For the layout of \c primary_variables_mat see collectDofsToMatrix().
template <typename NsAndWeight>
void interpolate(
    Eigen::MatrixXd const& primary_variables_mat,
    std::vector<std::reference_wrapper<ProcessVariable>> const& pv_refs,
    NsAndWeight const& ns_and_weight,
    Eigen::Ref<Eigen::VectorXd>
        interpolated_primary_variables)
{
    Eigen::Index component_flattened = 0;

    // We assume that all_process_variables_for_this_process have the same
    // order as the d.o.f. table. Therefore we can iterate over
    // all_process_variables_for_this_process.
    for (auto pv_ref : pv_refs)
    {
        auto const& pv = pv_ref.get();
        auto const num_comp = pv.getNumberOfGlobalComponents();
        auto const shp_fct_order = pv.getShapeFunctionOrder();
        auto const N = ns_and_weight.N(shp_fct_order);

        for (auto comp = decltype(num_comp){0}; comp < num_comp; ++comp)
        {
            // This computation assumes that there are no "holes" in the
            // primary_variables_mat. I.e., all nodal d.o.f. for a certain
            // (var, comp) must be stored contiguously in the respective column
            // of primary_variables_mat.
            interpolated_primary_variables[component_flattened] =
                N *
                primary_variables_mat.col(component_flattened).head(N.size());
            component_flattened++;
        }
    }
}
}  // namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
