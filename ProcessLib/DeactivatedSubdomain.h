/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DeactivatedSubdomain.h
 *
 * Created on November 29, 2018, 10:50 AM
 */
#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// TimeInterval cannot be forwardly declared because that
// std::unique_ptr<BaseLib::TimeInterval> type member requires its full
// definition (see https://stackoverflow.com/a/6089065).
#include "BaseLib/TimeInterval.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MeshLib
{
class Mesh;
class Node;
}  // namespace MeshLib

namespace ProcessLib
{
struct DeactivatedSubdomainMesh
{
    DeactivatedSubdomainMesh(
        std::unique_ptr<MeshLib::Mesh> deactivated_subdomain_mesh_,
        std::vector<MeshLib::Node*>&& inner_nodes_);

    std::unique_ptr<MeshLib::Mesh> const mesh;
    std::vector<MeshLib::Node*> const inner_nodes;
};

struct DeactivatedSubdomain
{
    DeactivatedSubdomain(
        std::unique_ptr<BaseLib::TimeInterval> time_interval_,
        std::vector<int>&& materialIDs_,
        std::vector<std::unique_ptr<DeactivatedSubdomainMesh>>&&
            deactivated_subdomain_meshes_);

    bool includesTimeOf(double const t) const;

    std::unique_ptr<BaseLib::TimeInterval const> const time_interval;

    /// The material IDs of the deactivated the subdomains
    std::vector<int> const materialIDs;

    std::vector<std::unique_ptr<DeactivatedSubdomainMesh>> const
        deactivated_subdomain_meshes;

    static const std::string zero_parameter_name;
};

std::vector<std::unique_ptr<DeactivatedSubdomain const>>
createDeactivatedSubdomains(BaseLib::ConfigTree const& config,
                            MeshLib::Mesh const& mesh);

}  // namespace ProcessLib
