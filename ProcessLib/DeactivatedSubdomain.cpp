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
    std::vector<std::size_t>&& bulk_element_ids_,
    std::vector<MeshLib::Node*>&& inner_nodes_,
    std::vector<MeshLib::Node*>&& outer_nodes_)
    : mesh(std::move(deactivated_subdomain_mesh_)),
      bulk_element_ids(std::move(bulk_element_ids_)),
      inner_nodes(std::move(inner_nodes_)),
      outer_nodes(std::move(outer_nodes_))
{
}

DeactivatedSubdomain::DeactivatedSubdomain(
    MathLib::PiecewiseLinearInterpolation time_interval_,
    std::optional<std::pair<Eigen::Vector3d, Eigen::Vector3d>>
        line_segment,
    std::unique_ptr<DeactivatedSubdomainMesh>&& deactivated_subdomain_mesh_,
    ParameterLib::Parameter<double> const* const boundary_value_parameter)
    : time_interval(std::move(time_interval_)),
      line_segment(line_segment),
      deactivated_subdomain_mesh(std::move(deactivated_subdomain_mesh_)),
      boundary_value_parameter(boundary_value_parameter)
{
}

bool DeactivatedSubdomain::isInTimeSupportInterval(double const t) const
{
    return time_interval.getSupportMin() <= t &&
           t <= time_interval.getSupportMax();
}

bool DeactivatedSubdomain::isDeactivated(MeshLib::Element const& element,
                                         double const time) const
{
    auto const& bulk_element_ids = deactivated_subdomain_mesh->bulk_element_ids;
    if (!BaseLib::contains(bulk_element_ids, element.getID()))
    {
        return false;
    }

    if (!line_segment)
    {
        return true;
    }

    auto const& element_center = getCenterOfGravity(element);
    // Line from a to b.
    auto const& a = line_segment->first;
    auto const& b = line_segment->second;
    // Tangent vector t = (b - a)/|b - a|.
    Eigen::Vector3d const t = (b - a).normalized();

    // Position r on the line at given time.
    auto const curve_position = time_interval.getValue(time);
    Eigen::Vector3d const r = a + t * curve_position;

    // Return true if p is "behind" the plane through r.
    return (element_center.asEigenVector3d() - r).dot(t) <= 0;
}
}  // namespace ProcessLib
