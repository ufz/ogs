/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DeactivatedSubdomain.h
 *
 * Created on November 29, 2018, 10:50 AM
 */
#pragma once

#include <Eigen/Dense>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MathLib/Point3d.h"

#include "processlib_export.h"

namespace MeshLib
{
class Element;
class Mesh;
class Node;
}  // namespace MeshLib

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib
{
struct DeactivatedSubdomainMesh
{
    DeactivatedSubdomainMesh(
        std::unique_ptr<MeshLib::Mesh> deactivated_subdomain_mesh_,
        std::vector<std::size_t>&& bulk_element_ids_,
        std::vector<MeshLib::Node*>&& inner_nodes_,
        std::vector<MeshLib::Node*>&& outer_nodes_);

    /// A mesh created from material ids (independent of time) for the
    /// deactivated subdomain.
    std::unique_ptr<MeshLib::Mesh> const mesh;
    std::vector<std::size_t> const bulk_element_ids;

    /// Inner nodes owned only by elements of the deactivated subdomain.
    /// \see ProcessLib::createDeactivatedSubdomainMesh()
    std::vector<MeshLib::Node*> const inner_nodes;
    /// Outer nodes owned by elements of the deactivated subdomain as well as
    /// other elements not being part of this deactivated subdomain mesh.
    /// \see ProcessLib::createDeactivatedSubdomainMesh()
    std::vector<MeshLib::Node*> const outer_nodes;
};

/// Time depend subdomain deactivation.
///
/// Subdomain deactivation is space and time dependent.
/// The spatial extent of deactivated elements is defined through a set of
/// material ids.
/// Additionally a line segment separates through its normal plane and position
/// on the line segment an active and inactive part.
///
/// The subdomain can be deactivated at once using a time interval.
/// For fine grained control a time curve can be specified.
/// It maps current time to position on the line segment.
/// Elements, which center points lie left of this position are deactivated.
///
/// \internal
/// The deactivated elements are excluded from the assembly, pre and post call,
/// secondary variables computation and the like. To keep the size of the global
/// linear equation system artificial Dirichlet boundary conditions are applied
/// on the interior of the deactivated subdomain. The nodes on the border
/// between the active and inactive elements are not affected.
struct DeactivatedSubdomain
{
    DeactivatedSubdomain(
        MathLib::PiecewiseLinearInterpolation time_interval_,
        std::optional<std::pair<Eigen::Vector3d, Eigen::Vector3d>>
            line_segment,
        std::unique_ptr<DeactivatedSubdomainMesh>&& deactivated_subdomain_mesh_,
        ParameterLib::Parameter<double> const* boundary_value_parameter);

    /// \returns true if the given time is included in the subdomains time
    /// support interval.
    bool isInTimeSupportInterval(double const t) const;

    /// \returns true if the element is in the deactivated part of the
    /// subdomain.
    /// If the line segment is available additionally the element's centre point
    /// is used to evaluate if the element is already, depending on time curve,
    /// active or inactive.  For this the domain is split into two parts by a
    /// plane defined as a normal plane of the line segment and the position on
    /// the line segment, where the latter is defined by the time curve.
    bool isDeactivated(MeshLib::Element const& element,
                       double const time) const;

    MathLib::PiecewiseLinearInterpolation const time_interval;

    /// Line segment along which excavation progresses. Represented by start and
    /// end points.
    std::optional<std::pair<Eigen::Vector3d, Eigen::Vector3d>> line_segment;

    std::unique_ptr<DeactivatedSubdomainMesh> const deactivated_subdomain_mesh;

    /// A pararameter for the optional Dirichlet boundary condition applied on
    /// the surface of the deactivated subdomain/excavation.
    ParameterLib::Parameter<double> const* const boundary_value_parameter;

    PROCESSLIB_EXPORT static const std::string zero_parameter_name;
};
}  // namespace ProcessLib
