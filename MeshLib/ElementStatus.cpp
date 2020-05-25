/**
 * \file
 * \author Karsten Rink
 * \date   2012-12-18
 * \brief  Implementation of the ElementStatus class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementStatus.h"

#include "Mesh.h"
#include "MeshLib/Node.h"
#include "Elements/Element.h"

namespace MeshLib {

ElementStatus::ElementStatus(Mesh const* const mesh, bool hasAnyInactive)
    : mesh_(mesh),
      element_status_(mesh->getNumberOfElements(), true),
      hasAnyInactive_(hasAnyInactive)
{
    const std::vector<MeshLib::Node*>& nodes(mesh_->getNodes());
    std::transform(
        begin(nodes), end(nodes), back_inserter(active_nodes_),
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
        for (auto material_id : vec_inactive_matIDs) {
            for (auto e : mesh_->getElements()) {
                if ((*materialIds)[e->getID()] == material_id) {
                    setElementStatus(e->getID(), false);
                }
            }
        }
    }

    vec_active_eles_.reserve(getNumberOfActiveElements());
    const std::size_t nElems (mesh_->getNumberOfElements());
    for (std::size_t i = 0; i < nElems; ++i)
    {
        if (element_status_[i])
        {
            vec_active_eles_.push_back(
                const_cast<MeshLib::Element*>(mesh_->getElement(i)));
        }
    }

    vec_active_nodes_.reserve(this->getNumberOfActiveNodes());
    const std::size_t nNodes (mesh_->getNumberOfNodes());
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        if (active_nodes_[i] > 0)
        {
            vec_active_nodes_.push_back(
                const_cast<MeshLib::Node*>(mesh_->getNode(i)));
        }
    }

    DBUG(
        "Deactivated {:d} materials and resulting active {:d} nodes and {:d} "
        "elements",
        vec_inactive_matIDs.size(),
        vec_active_nodes_.size(),
        vec_active_eles_.size());
}

std::vector<MeshLib::Element*> const& ElementStatus::getActiveElements() const
{
    if (hasAnyInactive_)
    {
        return vec_active_eles_;
    }

    return mesh_->getElements();
}

std::vector<MeshLib::Node*> const& ElementStatus::getActiveNodes() const
{
    if (hasAnyInactive_)
    {
        return vec_active_nodes_;
    }

    return mesh_->getNodes();
}

std::size_t ElementStatus::getNumberOfActiveNodes() const
{
    return active_nodes_.size() - std::count(active_nodes_.cbegin(), active_nodes_.cend(), 0);
}

std::size_t ElementStatus::getNumberOfActiveElements() const
{
    return static_cast<std::size_t>(
        std::count(element_status_.cbegin(), element_status_.cend(), true));
}

void ElementStatus::setElementStatus(std::size_t i, bool status)
{
    if (element_status_[i] != status)
    {
        const int change = (status) ? 1 : -1;
        element_status_[i] = status;
        const unsigned nElemNodes (mesh_->getElement(i)->getNumberOfNodes());
        MeshLib::Node const*const*const nodes = mesh_->getElement(i)->getNodes();
        for (unsigned j=0; j<nElemNodes; ++j)
        {
            assert(active_nodes_[j] < 255);  // if one node has >255 connected
                                             // elements the data type is too
                                             // small
            active_nodes_[nodes[j]->getID()] += change;
        }
    }
}

bool ElementStatus::isActiveNode(MeshLib::Node const* node) const
{
    return active_nodes_[node->getID()]>0;
}

}  // namespace MeshLib
