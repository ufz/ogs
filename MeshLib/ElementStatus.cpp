/**
 * \file
 * \author Karsten Rink
 * \date   2012-12-18
 * \brief  Implementation of the ElementStatus class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementStatus.h"

#include "Elements/Element.h"
#include "Mesh.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
ElementStatus::ElementStatus(Mesh const* const mesh, bool hasAnyInactive)
    : _mesh(mesh),
      _element_status(mesh->getNumberOfElements(), true),
      _hasAnyInactive(hasAnyInactive)
{
    const std::vector<MeshLib::Node*>& nodes(_mesh->getNodes());
    std::transform(
        begin(nodes), end(nodes), back_inserter(_active_nodes),
        [](Node const* const n) { return n->getNumberOfElements(); });
}

ElementStatus::ElementStatus(Mesh const* const mesh,
                             std::vector<int> const& vec_inactive_matIDs)
    : ElementStatus(mesh, !vec_inactive_matIDs.empty())
{
    if (mesh->getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        auto* const materialIds =
            mesh->getProperties().getPropertyVector<int>("MaterialIDs");
        for (auto material_id : vec_inactive_matIDs)
        {
            for (auto e : _mesh->getElements())
            {
                if ((*materialIds)[e->getID()] == material_id)
                {
                    setElementStatus(e->getID(), false);
                }
            }
        }
    }

    _vec_active_eles.reserve(getNumberOfActiveElements());
    const std::size_t nElems(_mesh->getNumberOfElements());
    for (std::size_t i = 0; i < nElems; ++i)
    {
        if (_element_status[i])
        {
            _vec_active_eles.push_back(
                const_cast<MeshLib::Element*>(_mesh->getElement(i)));
        }
    }

    _vec_active_nodes.reserve(this->getNumberOfActiveNodes());
    const std::size_t nNodes(_mesh->getNumberOfNodes());
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        if (_active_nodes[i] > 0)
        {
            _vec_active_nodes.push_back(
                const_cast<MeshLib::Node*>(_mesh->getNode(i)));
        }
    }

    DBUG(
        "Deactivated {:d} materials and resulting active {:d} nodes and {:d} "
        "elements",
        vec_inactive_matIDs.size(),
        _vec_active_nodes.size(),
        _vec_active_eles.size());
}

std::vector<MeshLib::Element*> const& ElementStatus::getActiveElements() const
{
    if (_hasAnyInactive)
    {
        return _vec_active_eles;
    }

    return _mesh->getElements();
}

std::vector<MeshLib::Node*> const& ElementStatus::getActiveNodes() const
{
    if (_hasAnyInactive)
    {
        return _vec_active_nodes;
    }

    return _mesh->getNodes();
}

std::size_t ElementStatus::getNumberOfActiveNodes() const
{
    return _active_nodes.size() -
           std::count(_active_nodes.cbegin(), _active_nodes.cend(), 0);
}

std::size_t ElementStatus::getNumberOfActiveElements() const
{
    return static_cast<std::size_t>(
        std::count(_element_status.cbegin(), _element_status.cend(), true));
}

void ElementStatus::setElementStatus(std::size_t i, bool status)
{
    if (_element_status[i] != status)
    {
        const int change = (status) ? 1 : -1;
        _element_status[i] = status;
        const unsigned nElemNodes(_mesh->getElement(i)->getNumberOfNodes());
        MeshLib::Node const* const* const nodes =
            _mesh->getElement(i)->getNodes();
        for (unsigned j = 0; j < nElemNodes; ++j)
        {
            assert(_active_nodes[j] < 255);  // if one node has >255 connected
                                             // elements the data type is too
                                             // small
            _active_nodes[nodes[j]->getID()] += change;
        }
    }
}

bool ElementStatus::isActiveNode(MeshLib::Node const* node) const
{
    return _active_nodes[node->getID()] > 0;
}

}  // namespace MeshLib
