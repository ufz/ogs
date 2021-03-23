/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace ProcessLib
{
const std::string DeactivatedSubdomain::zero_parameter_name =
    "zero_for_element_deactivation_approach";

DeactivatedSubdomainMesh::DeactivatedSubdomainMesh(
    std::unique_ptr<MeshLib::Mesh> deactivated_subdomain_mesh_,
    std::vector<MeshLib::Node*>&& inner_nodes_)
    : mesh(std::move(deactivated_subdomain_mesh_)),
      inner_nodes(std::move(inner_nodes_))
{
}

DeactivatedSubdomain::DeactivatedSubdomain(
    MathLib::PiecewiseLinearInterpolation time_interval_,
    std::pair<Eigen::Vector3d, Eigen::Vector3d>
        line_segment,
    std::vector<int>&& materialIDs_,
    std::vector<std::unique_ptr<DeactivatedSubdomainMesh>>&&
        deactivated_subdomain_meshes_)
    : time_interval(std::move(time_interval_)),
      line_segment(line_segment),
      materialIDs(std::move(materialIDs_)),
      deactivated_subdomain_meshes(std::move(deactivated_subdomain_meshes_))
{
}

bool DeactivatedSubdomain::includesTimeOf(double const t) const
{
    return time_interval.getSupportMin() <= t &&
           t <= time_interval.getSupportMax();
}

}  // namespace ProcessLib
