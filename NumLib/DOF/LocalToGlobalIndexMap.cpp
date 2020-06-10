/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalToGlobalIndexMap.h"

#include <algorithm>
#include <numeric>
#include <unordered_set>

namespace NumLib
{

namespace
{

// Make the cumulative sum of an array, which starts with zero
template <typename T>
std::vector<T> to_cumulative(std::vector<T> const& vec)
{
    std::vector<T> result(vec.size()+1, 0);
    std::partial_sum(vec.begin(), vec.end(), result.begin()+1);

    return result;
}

}  // namespace

int LocalToGlobalIndexMap::getGlobalComponent(int const variable_id,
                                              int const component_id) const
{
    return _variable_component_offsets[variable_id] + component_id;
}

template <typename ElementIterator>
void LocalToGlobalIndexMap::findGlobalIndicesWithElementID(
    ElementIterator first, ElementIterator last,
    std::vector<MeshLib::Node*> const& nodes, std::size_t const mesh_id,
    const int comp_id, const int comp_id_write)
{
    std::unordered_set<MeshLib::Node*> const set_nodes(nodes.begin(), nodes.end());

    // For each element find the global indices for node/element
    // components.
    for (ElementIterator e = first; e != last; ++e)
    {
        LineIndex indices;
        indices.reserve((*e)->getNumberOfNodes());

        for (auto* n = (*e)->getNodes();
             n < (*e)->getNodes()+(*e)->getNumberOfNodes(); ++n)
        {
            // Check if the element's node is in the given list of nodes.
            if (set_nodes.find(*n) == set_nodes.end())
            {
                continue;
            }
            MeshLib::Location l(
                mesh_id, MeshLib::MeshItemType::Node, (*n)->getID());
            indices.push_back(_mesh_component_map.getGlobalIndex(l, comp_id));
        }

        indices.shrink_to_fit();
        _rows((*e)->getID(), comp_id_write) = std::move(indices);
    }
}

template <typename ElementIterator>
void LocalToGlobalIndexMap::findGlobalIndices(
    ElementIterator first, ElementIterator last,
    std::vector<MeshLib::Node*> const& nodes, std::size_t const mesh_id,
    const int comp_id, const int comp_id_write)
{
    _rows.resize(std::distance(first, last), _mesh_subsets.size());

    std::unordered_set<MeshLib::Node*> const set_nodes(nodes.begin(), nodes.end());

    // For each element find the global indices for node/element
    // components.
    std::size_t elem_id = 0;
    for (ElementIterator e = first; e != last; ++e, ++elem_id)
    {
        LineIndex indices;
        indices.reserve((*e)->getNumberOfNodes());

        for (auto* n = (*e)->getNodes();
             n < (*e)->getNodes() + (*e)->getNumberOfNodes(); ++n)
        {
            // Check if the element's node is in the given list of nodes.
            if (set_nodes.find(*n) == set_nodes.end())
            {
                continue;
            }
            MeshLib::Location l(
                mesh_id, MeshLib::MeshItemType::Node, (*n)->getID());
            auto const global_index =
                _mesh_component_map.getGlobalIndex(l, comp_id);
            if (global_index == std::numeric_limits<GlobalIndexType>::max())
            {
                continue;
            }
            indices.push_back(global_index);
        }

        indices.shrink_to_fit();
        _rows(elem_id, comp_id_write) = std::move(indices);
    }
}

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubset>&& mesh_subsets,
    NumLib::ComponentOrder const order)
    : LocalToGlobalIndexMap(std::move(mesh_subsets),
                            std::vector<int>(mesh_subsets.size(), 1), order)
{
}

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubset>&& mesh_subsets,
    std::vector<int> const& vec_var_n_components,
    NumLib::ComponentOrder const order)
    : _mesh_subsets(std::move(mesh_subsets)),
      _mesh_component_map(_mesh_subsets, order),
      _variable_component_offsets(to_cumulative(vec_var_n_components))
{
    // For each element of that MeshSubset save a line of global indices.
    for (int variable_id = 0; variable_id < static_cast<int>(vec_var_n_components.size());
         ++variable_id)
    {
        for (int component_id = 0; component_id < static_cast<int>(vec_var_n_components[variable_id]);
             ++component_id)
        {
            auto const global_component_id =
                getGlobalComponent(variable_id, component_id);

            auto const& ms = _mesh_subsets[global_component_id];
            std::size_t const mesh_id = ms.getMeshID();

            findGlobalIndices(ms.elementsBegin(), ms.elementsEnd(),
                              ms.getNodes(), mesh_id, global_component_id,
                              global_component_id);
        }
    }
}

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubset>&& mesh_subsets,
    std::vector<int> const& vec_var_n_components,
    std::vector<std::vector<MeshLib::Element*> const*> const& vec_var_elements,
    NumLib::ComponentOrder const order)
    : _mesh_subsets(std::move(mesh_subsets)),
      _mesh_component_map(_mesh_subsets, order),
      _variable_component_offsets(to_cumulative(vec_var_n_components))
{
    assert(vec_var_n_components.size() == vec_var_elements.size());

    // For each element of that MeshSubset save a line of global indices.

    // _rows should be resized based on an element ID
    std::size_t max_elem_id = 0;
    for (std::vector<MeshLib::Element*>const* eles : vec_var_elements)
    {
        for (auto e : *eles)
        {
            max_elem_id = std::max(max_elem_id, e->getID());
        }
    }
    _rows.resize(max_elem_id + 1, _mesh_subsets.size());

    for (int variable_id = 0; variable_id < static_cast<int>(vec_var_n_components.size());
         ++variable_id)
    {
        std::vector<MeshLib::Element*> const& var_elements = *vec_var_elements[variable_id];
        for (int component_id = 0; component_id < static_cast<int>(vec_var_n_components[variable_id]);
             ++component_id)
        {
            auto const global_component_id =
                getGlobalComponent(variable_id, component_id);

            auto const& ms = _mesh_subsets[global_component_id];
            std::size_t const mesh_id = ms.getMeshID();

            findGlobalIndicesWithElementID(
                var_elements.cbegin(), var_elements.cend(), ms.getNodes(),
                mesh_id, global_component_id, global_component_id);
        }
    }
}

LocalToGlobalIndexMap::LocalToGlobalIndexMap(
    std::vector<MeshLib::MeshSubset>&& mesh_subsets,
    std::vector<int> const& global_component_ids,
    std::vector<int> const& variable_component_offsets,
    std::vector<MeshLib::Element*> const& elements,
    NumLib::MeshComponentMap&& mesh_component_map,
    LocalToGlobalIndexMap::ConstructorTag /*unused*/)
    : _mesh_subsets(std::move(mesh_subsets)),
      _mesh_component_map(std::move(mesh_component_map)),
      _variable_component_offsets(variable_component_offsets)
{
    // Each subset in the mesh_subsets represents a single component.
    if (_mesh_subsets.size() != global_component_ids.size())
    {
        OGS_FATAL(
            "Number of mesh subsets is not equal to number of components. "
            "There are {:d} mesh subsets and {:d} components.",
            _mesh_subsets.size(), global_component_ids.size());
    }

    for (int i = 0; i < static_cast<int>(global_component_ids.size()); ++i)
    {
        auto const& ms = _mesh_subsets[i];

        // For all MeshSubset in mesh_subsets and each element of that
        // MeshSubset save a line of global indices.
        std::size_t const mesh_id = ms.getMeshID();

        findGlobalIndices(elements.cbegin(), elements.cend(), ms.getNodes(),
                          mesh_id, global_component_ids[i], i);
    }
}

LocalToGlobalIndexMap* LocalToGlobalIndexMap::deriveBoundaryConstrainedMap(
    int const variable_id,
    std::vector<int> const& component_ids,
    MeshLib::MeshSubset&& new_mesh_subset) const
{
    DBUG("Construct reduced local to global index map.");

    if (component_ids.empty())
    {
        OGS_FATAL("Expected non-empty vector of component ids.");
    }

    // Elements of the new_mesh_subset's mesh.
    std::vector<MeshLib::Element*> const& elements =
        new_mesh_subset.getMesh().getElements();

    // Create a subset of the current mesh component map.
    std::vector<int> global_component_ids;

    transform(cbegin(component_ids), cend(component_ids),
              back_inserter(global_component_ids),
              [&](auto const component_id) {
                  return getGlobalComponent(variable_id, component_id);
              });

    auto mesh_component_map = _mesh_component_map.getSubset(
        _mesh_subsets, new_mesh_subset, global_component_ids);

    // Create copies of the new_mesh_subset for each of the global components.
    // The last component is moved after the for-loop.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    for (int i = 0; i < static_cast<int>(global_component_ids.size()) - 1; ++i)
    {
        all_mesh_subsets.emplace_back(new_mesh_subset);
    }
    all_mesh_subsets.emplace_back(std::move(new_mesh_subset));

    return new LocalToGlobalIndexMap(
        std::move(all_mesh_subsets), global_component_ids,
        _variable_component_offsets, elements, std::move(mesh_component_map),
        ConstructorTag{});
}

std::unique_ptr<LocalToGlobalIndexMap>
LocalToGlobalIndexMap::deriveBoundaryConstrainedMap(
    MeshLib::MeshSubset&& new_mesh_subset) const
{
    DBUG("Construct reduced local to global index map.");

    // Create a subset of the current mesh component map.
    std::vector<int> global_component_ids;

    for (int i = 0; i < getNumberOfComponents(); ++i)
    {
        global_component_ids.push_back(i);
    }

    // Elements of the new_mesh_subset's mesh.
    std::vector<MeshLib::Element*> const& elements =
        new_mesh_subset.getMesh().getElements();

    auto mesh_component_map = _mesh_component_map.getSubset(
        _mesh_subsets, new_mesh_subset, global_component_ids);

    // Create copies of the new_mesh_subset for each of the global components.
    // The last component is moved after the for-loop.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    for (int i = 0; i < static_cast<int>(global_component_ids.size()) - 1; ++i)
    {
        all_mesh_subsets.emplace_back(new_mesh_subset);
    }
    all_mesh_subsets.emplace_back(std::move(new_mesh_subset));

    return std::make_unique<LocalToGlobalIndexMap>(
        std::move(all_mesh_subsets), global_component_ids,
        _variable_component_offsets, elements, std::move(mesh_component_map),
        ConstructorTag{});
}

std::size_t LocalToGlobalIndexMap::dofSizeWithGhosts() const
{
    return _mesh_component_map.dofSizeWithGhosts();
}

std::size_t LocalToGlobalIndexMap::dofSizeWithoutGhosts() const
{
    return _mesh_component_map.dofSizeWithoutGhosts();
}

int LocalToGlobalIndexMap::getNumberOfVariables() const
{
    return static_cast<int>(_variable_component_offsets.size()) - 1;
}

int LocalToGlobalIndexMap::getNumberOfVariableComponents(int variable_id) const
{
    assert(variable_id < getNumberOfVariables());
    return _variable_component_offsets[variable_id + 1] -
           _variable_component_offsets[variable_id];
}

int LocalToGlobalIndexMap::getNumberOfComponents() const
{
    return _mesh_subsets.size();
}

std::size_t LocalToGlobalIndexMap::size() const
{
    return _rows.rows();
}

LocalToGlobalIndexMap::RowColumnIndices LocalToGlobalIndexMap::operator()(
    std::size_t const mesh_item_id, const int component_id) const
{
    return RowColumnIndices(_rows(mesh_item_id, component_id),
                            _columns(mesh_item_id, component_id));
}

std::size_t
LocalToGlobalIndexMap::getNumberOfElementDOF(std::size_t const mesh_item_id) const
{
    std::size_t ndof = 0;

    for (Table::Index c = 0; c < _rows.cols(); ++c)
    {
        ndof += _rows(mesh_item_id, c).size();
    }

    return ndof;
}

std::size_t
LocalToGlobalIndexMap::getNumberOfElementComponents(std::size_t const mesh_item_id) const
{
    std::size_t n = 0;
    for (Table::Index c = 0; c < _rows.cols(); ++c)
    {
        if (!_rows(mesh_item_id, c).empty())
        {
            n++;
        }
    }
    return n;
}

std::vector<int> LocalToGlobalIndexMap::getElementVariableIDs(
    std::size_t const mesh_item_id) const
{
    std::vector<int> vec;
    for (int i = 0; i < getNumberOfVariables(); i++)
    {
        for (int j=0; j<getNumberOfVariableComponents(i); j++)
        {
            auto comp_id = getGlobalComponent(i, j);
            if (!_rows(mesh_item_id, comp_id).empty())
            {
                vec.push_back(i);
            }
        }
    }
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());

    return vec;
}

GlobalIndexType LocalToGlobalIndexMap::getGlobalIndex(
    MeshLib::Location const& l,
    int const variable_id,
    int const component_id) const
{
    auto const c = getGlobalComponent(variable_id, component_id);
    return _mesh_component_map.getGlobalIndex(l, c);
}

GlobalIndexType LocalToGlobalIndexMap::getGlobalIndex(
    MeshLib::Location const& l, int const global_component_id) const
{
    return _mesh_component_map.getGlobalIndex(l, global_component_id);
}

/// Forwards the respective method from MeshComponentMap.
std::vector<GlobalIndexType> LocalToGlobalIndexMap::getGlobalIndices(
    const MeshLib::Location& l) const
{
    return _mesh_component_map.getGlobalIndices(l);
}

/// Get ghost indices, forwarded from MeshComponentMap.
std::vector<GlobalIndexType> const& LocalToGlobalIndexMap::getGhostIndices()
    const
{
    return _mesh_component_map.getGhostIndices();
}

/// Computes the index in a local (for DDC) vector for a given location and
/// component; forwarded from MeshComponentMap.
GlobalIndexType LocalToGlobalIndexMap::getLocalIndex(
    MeshLib::Location const& l, std::size_t const comp_id,
    std::size_t const range_begin, std::size_t const range_end) const
{
    return _mesh_component_map.getLocalIndex(l, comp_id, range_begin,
                                             range_end);
}

MeshLib::MeshSubset const& LocalToGlobalIndexMap::getMeshSubset(
    int const variable_id, int const component_id) const
{
    return getMeshSubset(getGlobalComponent(variable_id, component_id));
}

MeshLib::MeshSubset const& LocalToGlobalIndexMap::getMeshSubset(
    int const global_component_id) const
{
    return _mesh_subsets[global_component_id];
}

#ifndef NDEBUG
std::ostream& operator<<(std::ostream& os, LocalToGlobalIndexMap const& map)
{
    std::size_t const max_lines = 10;
    std::size_t lines_printed = 0;

    os << "Rows of the local to global index map; " << map._rows.size()
        << " rows\n";
    for (std::size_t e=0; e<map.size(); ++e)
    {
        os << "== e " << e << " ==\n";
        for (int c = 0; c < map.getNumberOfComponents(); ++c)
        {
            auto const& line = map._rows(e, c);

            os << "c" << c << " { ";
            std::copy(line.cbegin(), line.cend(),
                std::ostream_iterator<std::size_t>(os, " "));
            os << " }\n";
        }

        if (lines_printed++ > max_lines)
        {
            os << "...\n";
            break;
        }
    }

    os << "Mesh component map:\n" << map._mesh_component_map;
    return os;
}
#endif  // NDEBUG

}   // namespace NumLib
