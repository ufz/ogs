/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DeactivatedSubdomain.cpp
 *
 * Created on November 29, 2018, 10:50 AM
 */
#include "DeactivatedSubdomain.h"

#include <range/v3/algorithm/contains.hpp>

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
const std::string DeactivatedSubdomain::zero_parameter_name =
    "zero_for_element_deactivation_approach";

bool DeactivatedSubdomain::isInTimeSupportInterval(double const t) const
{
    return time_interval.getSupportMin() <= t &&
           t <= time_interval.getSupportMax();
}

bool DeactivatedSubdomain::isDeactivated(MeshLib::Element const& element,
                                         double const time) const
{
    auto const& bulk_element_ids = deactivated_subdomain_mesh.bulk_element_ids;
    if (!bulk_element_ids.contains(element.getID()))
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
