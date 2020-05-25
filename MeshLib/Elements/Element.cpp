/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Element class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Element.h"

#include "BaseLib/Logging.h"

#include "MathLib/GeometricBasics.h"
#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

#include "Line.h"

namespace MeshLib {

Element::Element(std::size_t id)
    : nodes_(nullptr), id_(id), content_(-1.0), neighbors_(nullptr)
{
}

Element::~Element()
{
    delete [] this->nodes_;
    delete [] this->neighbors_;
}

void Element::setNeighbor(Element* neighbor, unsigned const face_id)
{
    if (neighbor == this)
    {
        return;
    }

    this->neighbors_[face_id] = neighbor;
}

boost::optional<unsigned> Element::addNeighbor(Element* e)
{
    if (e == this || e == nullptr || e->getDimension() != this->getDimension())
    {
        return boost::optional<unsigned>();
    }

    if (this->hasNeighbor(e))
    {
        return boost::optional<unsigned>();
    }

    Node* face_nodes[3];
    const unsigned nNodes (this->getNumberOfBaseNodes());
    const unsigned eNodes (e->getNumberOfBaseNodes());
    const Node* const* e_nodes = e->getNodes();
    unsigned count(0);
    const unsigned dim (this->getDimension());
    for (unsigned i(0); i < nNodes; i++)
    {
        for (unsigned j(0); j < eNodes; j++)
        {
            if (nodes_[i] == e_nodes[j])
            {
                face_nodes[count] = nodes_[i];
                // increment shared nodes counter and check if enough nodes are similar to be sure e is a neighbour of this
                if ((++count)>=dim)
                {
                    neighbors_[ this->identifyFace(face_nodes) ] = e;
                    return boost::optional<unsigned>(e->identifyFace(face_nodes));
                }
            }
        }
    }

    return boost::optional<unsigned>();
}

MeshLib::Node Element::getCenterOfGravity() const
{
    const unsigned nNodes (this->getNumberOfBaseNodes());
    MeshLib::Node center(0,0,0);
    for (unsigned i=0; i<nNodes; ++i)
    {
        center[0] += (*nodes_[i])[0];
        center[1] += (*nodes_[i])[1];
        center[2] += (*nodes_[i])[2];
    }
    center[0] /= nNodes;
    center[1] /= nNodes;
    center[2] /= nNodes;
    return center;
}

void Element::computeSqrEdgeLengthRange(double &min, double &max) const
{
    min = std::numeric_limits<double>::max();
    max = 0;
    const unsigned nEdges (this->getNumberOfEdges());
    for (unsigned i=0; i<nEdges; i++)
    {
        const double dist (MathLib::sqrDist(*getEdgeNode(i,0), *getEdgeNode(i,1)));
        min = (dist<min) ? dist : min;
        max = (dist>max) ? dist : max;
    }
}

void Element::computeSqrNodeDistanceRange(double &min, double &max, bool check_allnodes) const
{
    min = std::numeric_limits<double>::max();
    max = 0;
    const unsigned nnodes = check_allnodes ? getNumberOfNodes() : getNumberOfBaseNodes();
    for (unsigned i=0; i<nnodes; i++)
    {
        for (unsigned j=i+1; j<nnodes; j++)
        {
            const double dist (MathLib::sqrDist(*getNode(i), *getNode(j)));
            min = std::min(dist, min);
            max = std::max(dist, max);
        }
    }
}

const Element* Element::getNeighbor(unsigned i) const
{
#ifndef NDEBUG
    if (i < getNumberOfNeighbors())
#endif
    {
        return neighbors_[i];
    }
#ifndef NDEBUG
    ERR("Error in MeshLib::Element::getNeighbor() - Index does not exist.");
    return nullptr;
#endif
}

unsigned Element::getNodeIDinElement(const MeshLib::Node* node) const
{
    const unsigned nNodes (this->getNumberOfNodes());
    for (unsigned i(0); i < nNodes; i++)
    {
        if (node == nodes_[i])
        {
            return i;
        }
    }
    return std::numeric_limits<unsigned>::max();
}

const Node* Element::getNode(unsigned i) const
{
#ifndef NDEBUG
    if (i < getNumberOfNodes())
#endif
    {
        return nodes_[i];
    }
#ifndef NDEBUG
    ERR("Error in MeshLib::Element::getNode() - Index {:d} in {:s}", i,
        MeshElemType2String(getGeomType()));
    return nullptr;
#endif
}

void Element::setNode(unsigned idx, Node* node)
{
#ifndef NDEBUG
    if (idx < getNumberOfNodes())
#endif
    {
        nodes_[idx] = node;
    }
}

std::size_t Element::getNodeIndex(unsigned i) const
{
#ifndef NDEBUG
    if (i < getNumberOfNodes())
#endif
    {
        return nodes_[i]->getID();
    }
#ifndef NDEBUG
    ERR("Error in MeshLib::Element::getNodeIndex() - Index does not exist.");
    return std::numeric_limits<std::size_t>::max();
#endif
}

bool Element::hasNeighbor(Element* elem) const
{
    unsigned nNeighbors (this->getNumberOfNeighbors());
    for (unsigned i = 0; i < nNeighbors; i++)
    {
        if (this->neighbors_[i] == elem)
        {
            return true;
        }
    }
    return false;
}

bool Element::isBoundaryElement() const
{
    return std::any_of(neighbors_, neighbors_ + this->getNumberOfNeighbors(),
        [](MeshLib::Element const*const e){ return e == nullptr; });
}

#ifndef NDEBUG
std::ostream& operator<<(std::ostream& os, Element const& e)
{
    os << "Element #" << e.id_ << " @ " << &e << " with " << e.getNumberOfNeighbors()
       << " neighbours\n";

    unsigned const nnodes = e.getNumberOfNodes();
    MeshLib::Node* const* const nodes = e.getNodes();
    os << "MeshElemType: "
        << static_cast<std::underlying_type<MeshElemType>::type>(e.getGeomType())
        << " with " << nnodes << " nodes: { ";
    for (unsigned n = 0; n < nnodes; ++n)
    {
        os << nodes[n]->getID() << " @ " << nodes[n] << "  ";
    }
    os << "}\n";
    return os;
}
#endif  // NDEBUG

bool isPointInElementXY(MathLib::Point3d const& p, Element const& e)
{
    for(std::size_t i(0); i<e.getNumberOfBaseNodes(); ++i) {
        if (MathLib::sqrDist2d(p, *e.getNode(i)) <
            std::numeric_limits<double>::epsilon())
        {
            return true;
        }
    }

    if (e.getGeomType() == MeshElemType::TRIANGLE)
    {
        MathLib::Point3d const& n0(*e.getNode(0));
        MathLib::Point3d const& n1(*e.getNode(1));
        MathLib::Point3d const& n2(*e.getNode(2));

        return MathLib::isPointInTriangleXY(p, n0, n1, n2);
    }
    if (e.getGeomType() == MeshElemType::QUAD)
    {
        MathLib::Point3d const& n0(*e.getNode(0));
        MathLib::Point3d const& n1(*e.getNode(1));
        MathLib::Point3d const& n2(*e.getNode(2));
        MathLib::Point3d const& n3(*e.getNode(3));

        return MathLib::isPointInTriangleXY(p, n0, n1, n2) ||
               MathLib::isPointInTriangleXY(p, n0, n2, n3);
    }

    WARN("isPointInElementXY: element type '{:s}' is not supported.",
         MeshLib::MeshElemType2String(e.getGeomType()));
    return false;
}
}  // namespace MeshLib
