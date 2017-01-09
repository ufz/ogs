/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"
#include "MeshSubset.h"

namespace MeshLib
{

/// Collection of mesh subsets.
class MeshSubsets
{
public:

    /// Single mesh subset constructor.
    MeshSubsets(const MeshSubset* mesh_subset)
    {
        _mesh_subsets.push_back(mesh_subset);
        _n_total_items = mesh_subset->getNumberOfTotalItems();
    }

    /// Construct MeshSubsets from a range of MeshSubset. InputIterator must
    /// dereference to MeshSubset*.
    /// \pre All meshes of each of the MeshSubset objects must unique,
    /// an exception is thrown otherwise.
    template <typename InputIterator>
    MeshSubsets(InputIterator const& first, InputIterator const& last)
        : _mesh_subsets(first, last)
    {
        if (!areMeshSubsetMeshesUnique())
        {
            OGS_FATAL("Mesh ids of input mesh subsets are not unique.");
        }

        _n_total_items = std::accumulate(first, last, 0u,
            [](std::size_t const& sum, MeshSubset const* const mesh_subset)
            {
                return sum + mesh_subset->getNumberOfTotalItems();
            });
    }

    /// return the total number of mesh items (in all meshes) where this component is assigned
    std::size_t getNumberOfMeshItems() const
    {
        return _n_total_items;
    }

    /// Number of saved mesh subsets.
    std::size_t size() const
    {
        return _mesh_subsets.size();
    }

    /// return MeshSubset
    const MeshSubset& getMeshSubset(std::size_t mesh_index) const
    {
        return *_mesh_subsets[mesh_index];
    }

    std::vector<const MeshSubset*>::const_iterator begin() const
    {
        return _mesh_subsets.begin();
    }

    std::vector<const MeshSubset*>::const_iterator end() const
    {
        return _mesh_subsets.end();
    }
private:
    /// returns true if all mesh ids of _mesh_subsets elements are different.
    bool areMeshSubsetMeshesUnique() const
    {
        std::vector<std::size_t> mesh_ids;
        std::transform(_mesh_subsets.cbegin(), _mesh_subsets.cend(),
            std::back_inserter(mesh_ids), std::mem_fun(&MeshSubset::getMeshID));

        std::sort(mesh_ids.begin(), mesh_ids.end());
        auto const it = std::adjacent_find(mesh_ids.cbegin(), mesh_ids.cend());

        return mesh_ids.cend() == it;
    }

private:
    std::vector<const MeshSubset*> _mesh_subsets;

    /// Number of all mesh entities on all subsets.
    std::size_t _n_total_items;
};

}   // namespace MeshLib
