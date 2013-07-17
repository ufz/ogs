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

#ifndef MESHSUBSETS_H_
#define MESHSUBSETS_H_

#include <vector>
#include <numeric>

#include "MeshSubset.h"

namespace VecMatOnMeshLib
{

/// Collection of mesh subsets.
class MeshSubsets
{
public:

    /// Single mesh subset constructor.
    explicit MeshSubsets(const MeshSubset* mesh_subset)
    {
        _mesh_subsets.push_back(mesh_subset);
        _n_total_items = mesh_subset->getNTotalItems();
    }

	/// Construct MeshSubsets from a range of MeshSubset. InputIterator must
	/// dereferece to MeshSubset*.
	template <typename InputIterator>
    explicit MeshSubsets(InputIterator const& first, InputIterator const& last)
		: _mesh_subsets(first, last)
    {
        _n_total_items = std::accumulate(first, last, 0u,
			[](std::size_t const& sum, MeshSubset const* const mesh_subset)
			{
				return sum + mesh_subset->getNTotalItems();
			});
    }

    /// return the total number of mesh items (in all meshes) where this component is assigned
    std::size_t getNMeshItems() const
	{
		return _n_total_items;
	}

    /// return the number of related meshes
    unsigned getNMeshes() const
	{
		return _mesh_subsets.size();
	}

    /// return MeshSubset
    const MeshSubset& getMeshItems(std::size_t mesh_index) const
	{
		return *_mesh_subsets[mesh_index];
	}

private:
    std::vector<const MeshSubset*> _mesh_subsets;

	/// Number of all mesh entities on all subsets.
    std::size_t _n_total_items;
};

}

#endif	// MESHSUBSETS_H_
