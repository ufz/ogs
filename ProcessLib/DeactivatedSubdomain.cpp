/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DeactivatedSubdomain.cpp
 *
 * Created on November 29, 2018, 10:50 AM
 */
#include "DeactivatedSubdomain.h"

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
const std::string DeactivatedSubdomain::zero_parameter_name =
    "zero_for_element_deactivation_approach";

DeactivatedSubdomainMesh::DeactivatedSubdomainMesh(
    std::unique_ptr<MeshLib::Mesh> deactivated_subdomain_mesh_,
    std::vector<MeshLib::Node*>&& inner_nodes_,
    std::vector<MeshLib::Node*>&& outer_nodes_)
    : mesh(std::move(deactivated_subdomain_mesh_)),
      inner_nodes(std::move(inner_nodes_)),
      outer_nodes(std::move(outer_nodes_))
{
}

DeactivatedSubdomain::DeactivatedSubdomain(
    MathLib::PiecewiseLinearInterpolation time_interval_,
    std::pair<Eigen::Vector3d, Eigen::Vector3d>
        line_segment,
    std::vector<int>&& materialIDs_,
    std::vector<std::unique_ptr<DeactivatedSubdomainMesh>>&&
        deactivated_subdomain_meshes_,
    ParameterLib::Parameter<double> const* const boundary_value_parameter)
    : time_interval(std::move(time_interval_)),
      line_segment(line_segment),
      materialIDs(std::move(materialIDs_)),
      deactivated_subdomain_meshes(std::move(deactivated_subdomain_meshes_)),
      boundary_value_parameter(boundary_value_parameter)
{
}

bool DeactivatedSubdomain::isInTimeSupportInterval(double const t) const
{
    return time_interval.getSupportMin() <= t &&
           t <= time_interval.getSupportMax();
}

bool DeactivatedSubdomain::isDeactivated(MathLib::Point3d const& point,
                                         double const time) const
{
    // Line from a to b.
    auto const& a = line_segment.first;
    auto const& b = line_segment.second;
    // Tangent vector t = (b - a)/|b - a|.
    Eigen::Vector3d const t = (b - a).normalized();

    // Position r on the line at given time.
    Eigen::Vector3d const r = a + t * time_interval.getValue(time);
    Eigen::Map<Eigen::Vector3d const> const p{point.data(), 3};

    // Return true if p is "behind" the plane through r.
    return (p - r).dot(t) <= 0;
}
}  // namespace ProcessLib
