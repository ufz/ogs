/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef COMPONENTDISTRIBUTION_H_
#define COMPONENTDISTRIBUTION_H_

#include <vector>
#include <numeric>

#include "MeshSubset.h"

namespace VecMatOnMeshLib
{

/**
 * Distribution information of a single data component
 *
 * This class contains information about on which mesh items a data component is
 * assigned.
 */
class ComponentDistribution
{
public:

    /**
     * constructor for a single mesh use
     *
     * This data component is distributed over the given mesh items in a single mesh
     * @param mesh_items
     */
    explicit ComponentDistribution(const MeshSubset* mesh_items)
    {
        _mesh_items.push_back(mesh_items);
        _n_total_items = mesh_items->getNTotalItems();
    }

    /**
     * constructor for multiple-mesh use
     *
     * This data component is distributed over the given mesh items in multiple meshes
     * @param vec_mesh_items   a vector of MeshSubset
     */
    explicit ComponentDistribution(const std::vector<MeshSubset*> &vec_mesh_items)
    : _mesh_items(vec_mesh_items.begin(), vec_mesh_items.end())
    {
        _n_total_items = std::accumulate(vec_mesh_items.begin(), vec_mesh_items.end(),
                                            0u,
                                            [](std::size_t sum, const MeshSubset* items)
                                            {
                                                return sum+items->getNTotalItems();
                                            }
                                           );
    }

    /// return the total number of mesh items (in all meshes) where this component is assigned
    std::size_t getNMeshItems() const { return _n_total_items; }

    /// return the number of related meshes
    unsigned getNMeshes() const { return _mesh_items.size(); }

    /// return MeshSubset
    const MeshSubset& getMeshItems(std::size_t mesh_index) const { return *_mesh_items[mesh_index]; }

private:
    std::vector<const MeshSubset*> _mesh_items;
    std::size_t _n_total_items;
};

}

#endif /* COMPONENTDISTRIBUTION_H_ */
