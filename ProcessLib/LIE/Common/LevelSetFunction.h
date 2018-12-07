/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <unordered_map>
#include <vector>

#include <Eigen/Eigen>

namespace ProcessLib
{
namespace LIE
{
struct FractureProperty;
struct JunctionProperty;

/// calculate the enrichment function for displacements at the given point
/// Remarks:
/// * branch/junction intersections of two fractures are supported in 2D
///
/// @parm frac_props fracture properties
/// @parm junction_props junction properties
/// @parm fracID_to_local a mapping table from a fracture ID to a local index in
/// frac_props
/// @parm x evaluating point coordiates
/// @return a vector of enrichment values for displacements
std::vector<double> u_global_enrichments(
    std::vector<FractureProperty*> const& frac_props,
    std::vector<JunctionProperty*> const& junction_props,
    std::unordered_map<int, int> const& fracID_to_local,
    Eigen::Vector3d const& x);

/// calculate the enrichment function for fracture relative displacements
/// Remarks:
/// * branch/junction intersections of two fractures are supported in 2D
///
/// @parm this_fracID the fracture ID
/// @parm frac_props fracture properties
/// @parm junction_props junction properties
/// @parm fracID_to_local a mapping table from a fracture ID to a local index in
/// frac_props
/// @parm x evaluating point coordiates
/// @return a vector of enrichment values for fracture relative displacements
std::vector<double> du_global_enrichments(
    std::size_t this_fracID,
    std::vector<FractureProperty*> const& frac_props,
    std::vector<JunctionProperty*> const& junction_props,
    std::unordered_map<int, int> const& fracID_to_local,
    Eigen::Vector3d const& x);

}  // namespace LIE
}  // namespace ProcessLib
