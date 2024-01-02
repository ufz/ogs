/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementSearch.h"

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/take.hpp>

#include "BaseLib/Algorithm.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
ElementSearch::ElementSearch(const MeshLib::Mesh& mesh) : _mesh(mesh) {}

template <typename Container, typename Predicate>
std::vector<std::size_t> filter(Container const& container, Predicate const& p)
{
    return ranges::views::filter(container, p) | views::ids |
           ranges::to<std::vector>;
}

std::size_t ElementSearch::searchByElementType(MeshElemType eleType)
{
    auto matchedIDs = filter(_mesh.getElements(), [&](MeshLib::Element const* e)
                             { return e->getGeomType() == eleType; });

    this->updateUnion(matchedIDs);
    return matchedIDs.size();
}

std::size_t ElementSearch::searchByContent(double eps)
{
    auto matchedIDs =
        filter(_mesh.getElements(), [&eps](MeshLib::Element const* e)
               { return e->getContent() < eps; });

    this->updateUnion(matchedIDs);
    return matchedIDs.size();
}

std::size_t ElementSearch::searchByBoundingBox(GeoLib::AABB const& aabb,
                                               bool const invert)
{
    auto matchedIDs = filter(
        _mesh.getElements(),
        [&aabb, invert](MeshLib::Element const* e)
        {
            // any node of element is in aabb.
            return ranges::any_of(
                e->nodes() | ranges::views::take(e->getNumberOfBaseNodes()),
                [&aabb, invert](auto const* n)
                { return (aabb.containsPoint(*n, 0) != invert); });
        });

    this->updateUnion(matchedIDs);
    return matchedIDs.size();
}

std::size_t ElementSearch::searchByNodeIDs(
    const std::vector<std::size_t>& nodes)
{
    std::vector<std::size_t> connected_elements;
    for (std::size_t node_id : nodes)
    {
        auto const& elements = _mesh.getElementsConnectedToNode(node_id);
        std::transform(begin(elements), end(elements),
                       back_inserter(connected_elements),
                       [](Element const* const e) { return e->getID(); });
    }

    BaseLib::makeVectorUnique(connected_elements);

    this->updateUnion(connected_elements);
    return connected_elements.size();
}

void ElementSearch::updateUnion(const std::vector<std::size_t>& vec)
{
    std::vector<std::size_t> vec_temp(vec.size() + _marked_elements.size());
    auto it = std::set_union(vec.begin(), vec.end(), _marked_elements.begin(),
                             _marked_elements.end(), vec_temp.begin());
    vec_temp.resize(it - vec_temp.begin());
    _marked_elements.assign(vec_temp.begin(), vec_temp.end());
}

}  // end namespace MeshLib
