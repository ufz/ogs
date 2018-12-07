/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

namespace ProcessLib
{
namespace LIE
{
struct JunctionProperty
{
    int junction_id;
    int node_id;
    Eigen::Vector3d coords;
    std::array<int, 2> fracture_IDs;

    virtual ~JunctionProperty() = default;
};

inline void setJunctionProperty(int junction_id,
                                MeshLib::Node const& junctionNode,
                                std::vector<int>
                                    frac_ids,
                                JunctionProperty& junction)
{
    junction.junction_id = junction_id;
    junction.node_id = junctionNode.getID();
    junction.coords = Eigen::Vector3d(junctionNode.getCoords());
    for (int j = 0; j < 2; j++)
        junction.fracture_IDs[j] = frac_ids[j];
}

}  // namespace LIE
}  // namespace ProcessLib
