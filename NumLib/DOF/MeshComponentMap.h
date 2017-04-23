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

#include "ComponentGlobalIndexDict.h"

namespace MeshLib
{
    class MeshSubsets;
}

namespace NumLib
{

/// Ordering of components in global matrix/vector.
enum class ComponentOrder
{
    BY_COMPONENT,   ///< Ordering data by component type
    BY_LOCATION     ///< Ordering data by spatial location
};

/// Multidirectional mapping between mesh entities and degrees of freedom.
class MeshComponentMap final
{
public:
    typedef MeshLib::Location Location;
public:
    /// \param components   a vector of components
    /// \param order        type of ordering values in a vector
    MeshComponentMap(
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>> const& components,
        ComponentOrder order);

    /// Creates a multi-component subset of the current mesh component map.
    /// The order (BY_LOCATION/BY_COMPONENT) of components is the same as of the
    /// current map.
    ///
    /// \note For each component the same mesh subset will be used.
    ///
    /// \param component_ids  The vector of global components id.
    /// \param components components that should remain in the created subset
    MeshComponentMap getSubset(std::vector<std::size_t> const& component_ids,
                               MeshLib::MeshSubsets const& components) const;

    /// The number of dofs including the those located in the ghost nodes.
    std::size_t dofSizeWithGhosts() const
    {
        return _dict.size();
    }

    /// Component ids at given location \c l.
    ///
    /// | Location | ComponentID |
    /// | -------- | ----------- |
    /// | l        | comp_id_1   |
    /// | l        |  ...        |
    /// | l        | comp_id_n   |
    std::vector<std::size_t> getComponentIDs(const Location &l) const;

    /// Global index of the given component id at given location \c l.
    ///
    /// | Location | ComponentID | GlobalIndex |
    /// | -------- | ----------- | ----------- |
    /// | l        | comp_id     | gi          |
    GlobalIndexType getGlobalIndex(Location const &l, std::size_t const comp_id) const;

    /// Global indices for all components at the given location \c l.
    ///
    /// If there is more than one component at the given location, the
    /// function returns a vector containing global indices for each component.
    ///
    /// | Location | ComponentID | GlobalIndex |
    /// | -------- | ----------- | ----------- |
    /// | l        | comp_id_1   | gi23        |
    /// | ...      |  ...        | ...         |
    /// | l        | comp_id_k   | gi45        |
    std::vector<GlobalIndexType> getGlobalIndices(const Location &l) const;

    /// Global indices for all components at all given locations \c ls ordered
    /// by location. The return list is sorted first by location.
    ///
    /// | Location | ComponentID | GlobalIndex |
    /// | -------- | ----------- | ----------- |
    /// | l_1      | comp_id_1   | gi23        |
    /// | ...      |  ...        | ...         |
    /// | l_1      | comp_id_k   | gi45        |
    /// | l_2      | comp_id_1   | gi46        |
    /// | ...      |  ...        | ...         |
    /// | l_2      | comp_id_m   | gi67        |
    /// | ...      |  ...        | ...         |
    /// | l_n      | comp_id_n   | gi78        |
    std::vector<GlobalIndexType> getGlobalIndicesByLocation(
        const std::vector<Location>& ls) const;

    /// Global indices for all components at all given locations \c ls ordered
    /// by component ids. The return list is sorted first by component ids.
    /// | Location | ComponentID | GlobalIndex |
    /// | -------- | ----------- | ----------- |
    /// | l_1      | comp_id_1   | gi23        |
    /// | ...      |  ...        | ...         |
    /// | l_k      | comp_id_1   | gi45        |
    /// | l_1      | comp_id_2   | gi46        |
    /// | ...      |  ...        | ...         |
    /// | l_m      | comp_id_2   | gi78        |
    /// | ...      |  ...        | ...         |
    /// | l_n      | comp_id_n   | gi89        |
    std::vector<GlobalIndexType> getGlobalIndicesByComponent(
        const std::vector<Location>& ls) const;

    /// Get the number of local unknowns excluding those associated
    /// with ghost nodes (for DDC with node-wise mesh partitioning).
    std::size_t dofSizeWithoutGhosts() const
    {
        return _num_local_dof;
    }

    /// Get ghost indices (for DDC).
    std::vector<GlobalIndexType> const& getGhostIndices() const
    {
        return _ghosts_indices;
    }

    /// Computes the index in a local (for DDC) vector for a given location and
    /// component.
    /// When domain decomposition is not used, it is equal to getGlobalIndex().
    /// The range is needed to compute the offset for non-ghost locations and
    /// also to map ghost locations.
    GlobalIndexType getLocalIndex(Location const& l, std::size_t const comp_id,
                                  std::size_t const range_begin,
                                  std::size_t const range_end) const;

    /// A value returned if no global index was found for the requested
    /// location/component. The value is implementation dependent.
    static GlobalIndexType const nop;

#ifndef NDEBUG
    const detail::ComponentGlobalIndexDict& getDictionary() const
    {
        return _dict;
    }

    friend std::ostream& operator<<(std::ostream& os, MeshComponentMap const& m)
    {
        os << "Dictionary size: " << m._dict.size() << "\n";
        for (auto l : m._dict)
            os << l << "\n";
        return os;
    }
#endif  // NDEBUG

private:
    /// Private constructor used by internally created mesh component maps.
    MeshComponentMap(detail::ComponentGlobalIndexDict& dict,
                     unsigned const num_components)
        : _dict(dict), _num_components(num_components)
    { }

    /// Looks up if a line is already stored in the dictionary.
    /// \attention The line for the location l and component id must exist,
    /// the behaviour is undefined otherwise.
    /// \return a copy of the line.
    detail::Line getLine(Location const& l, std::size_t const component_id) const;

    void renumberByLocation(GlobalIndexType offset=0);

    detail::ComponentGlobalIndexDict _dict;

    /// Number of local unknowns excluding those associated
    /// with ghost nodes (for domain decomposition).
    std::size_t _num_local_dof  = 0;

#ifdef USE_PETSC
    /// Number of global unknowns. Used internally only.
    std::size_t _num_global_dof = 0;
#endif

    /// Number of components
    /// introduced mainly for error checking
    unsigned const _num_components;

    /// Global ID for ghost entries
    std::vector<GlobalIndexType> _ghosts_indices;
};

}   // namespace NumLib
