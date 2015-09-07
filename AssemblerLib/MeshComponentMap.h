/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_MESHCOMPONENTMAP_H_
#define ASSEMBLERLIB_MESHCOMPONENTMAP_H_

#include <vector>

#include "MeshLib/Location.h"

#include "ComponentGlobalIndexDict.h"

namespace MeshLib
{
    class MeshSubsets;
}

namespace AssemblerLib
{

/// Ordering of components in global matrix/vector.
enum class ComponentOrder
{
    BY_COMPONENT,   ///< Ordering data by component type
    BY_LOCATION     ///< Ordering data by spatial location
};

/// Multidirectional mapping between mesh entities and degrees of freedom.
class MeshComponentMap
{
public:
    typedef MeshLib::Location Location;
public:
    /// \param components   a vector of components
    /// \param order        type of ordering values in a vector
    MeshComponentMap(std::vector<MeshLib::MeshSubsets*> const& components,
        ComponentOrder order);

    /// Creates a subset of the current mesh component map.
    /// The order of components is the same as of the current map.
    MeshComponentMap getSubset(
        std::vector<MeshLib::MeshSubsets*> const& components) const;

    /// The number of components in the map.
    std::size_t size() const
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
    std::size_t getGlobalIndex(Location const &l, std::size_t const comp_id) const;

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
    std::vector<std::size_t> getGlobalIndices(const Location &l) const;

    /// Global indices for all components at all given locations \c ls ordered
    /// as required by the template parameter ORDER.
    ///
    /// In case ORDER is ComponentOrder::BY_LOCATION the return list is sorted
    /// first by location.
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
    ///
    /// In case ORDER is ComponentOrder::BY_COMPONENT the return list is sorted
    /// first by component ids.
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
    template <ComponentOrder ORDER>
    std::vector<std::size_t> getGlobalIndices(const std::vector<Location> &ls) const;

    /// Get global indices for the given component
    ///
    /// Given the global indices for all compontents of a current mesh element,
    /// and the index of a component, this method will return the global indices
    /// of the component.
    ///
    /// \param cnt collection of global indices, e.g. for a mesh element
    /// \param component_id the component of interest
    std::vector<std::size_t>
    getIndicesForComponent(const std::vector<std::size_t>& cnt,
                           const unsigned component_id) const;


    unsigned getNumComponents() const { return _num_components; }

    /// A value returned if no global index was found for the requested
    /// location/component. The value is implementation dependent.
    static std::size_t const nop;

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
                     ComponentOrder const order,
                     unsigned const num_components)
        : _dict(dict), _order(order), _num_components(num_components)
    { }

    /// Looks up if a line is already stored in the dictionary.
    /// \attention The line for the location l and component id must exist,
    /// the behaviour is undefined otherwise.
    /// \return a copy of the line.
    detail::Line getLine(Location const& l, std::size_t const component_id) const;

    void renumberByLocation(std::size_t offset=0);

private:
    detail::ComponentGlobalIndexDict _dict;
    ComponentOrder const _order;
    unsigned const _num_components;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_MESHCOMPONENTMAP_H_
