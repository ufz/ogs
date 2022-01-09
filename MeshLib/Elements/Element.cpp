/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Element class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Element.h"

#include "BaseLib/Logging.h"
#include "Line.h"
#include "MathLib/GeometricBasics.h"
#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
Element::Element(std::size_t id) : _id(id), _neighbors(nullptr) {}

Element::~Element()
{
    delete[] this->_neighbors;
}

void Element::setNeighbor(Element* neighbor, unsigned const face_id)
{
    if (neighbor == this)
    {
        return;
    }

    this->_neighbors[face_id] = neighbor;
}

std::optional<unsigned> Element::addNeighbor(Element* e)
{
    const unsigned dim(this->getDimension());
    if (e == this || e == nullptr || e->getDimension() != dim)
    {
        return std::optional<unsigned>();
    }

    if (areNeighbors(this, e))
    {
        return std::optional<unsigned>();
    }

    Node const* face_nodes[3];
    const unsigned nNodes(this->getNumberOfBaseNodes());
    const unsigned eNodes(e->getNumberOfBaseNodes());
    const Node* const* e_nodes = e->getNodes();
    unsigned count(0);
    for (unsigned i(0); i < nNodes; i++)
    {
        for (unsigned j(0); j < eNodes; j++)
        {
            if (getNode(i) == e_nodes[j])
            {
                face_nodes[count] = getNode(i);
                // increment shared nodes counter and check if enough nodes are
                // similar to be sure e is a neighbour of this
                if ((++count) >= dim)
                {
                    _neighbors[this->identifyFace(face_nodes)] = e;
                    return std::optional<unsigned>(e->identifyFace(face_nodes));
                }
            }
        }
    }

    return std::optional<unsigned>();
}

bool Element::isBoundaryElement() const
{
    return std::any_of(_neighbors, _neighbors + this->getNumberOfNeighbors(),
                       [](MeshLib::Element const* const e)
                       { return e == nullptr; });
}

#ifndef NDEBUG
std::ostream& operator<<(std::ostream& os, Element const& e)
{
    os << "Element #" << e._id << " @ " << &e << " with "
       << e.getNumberOfNeighbors() << " neighbours\n";

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

bool areNeighbors(Element const* const element, Element const* const other)
{
    unsigned nNeighbors(element->getNumberOfNeighbors());
    for (unsigned i = 0; i < nNeighbors; i++)
    {
        if (element->getNeighbor(i) == other)
        {
            return true;
        }
    }
    return false;
}

bool hasZeroVolume(MeshLib::Element const& element)
{
    return element.getContent() < std::numeric_limits<double>::epsilon();
}

MeshLib::Node getCenterOfGravity(Element const& element)
{
    const unsigned nNodes(element.getNumberOfBaseNodes());
    MeshLib::Node center(0, 0, 0);
    for (unsigned i = 0; i < nNodes; ++i)
    {
        center[0] += (*element.getNode(i))[0];
        center[1] += (*element.getNode(i))[1];
        center[2] += (*element.getNode(i))[2];
    }
    center[0] /= nNodes;
    center[1] /= nNodes;
    center[2] /= nNodes;
    return center;
}

std::pair<double, double> computeSqrNodeDistanceRange(
    MeshLib::Element const& element, bool const check_allnodes)
{
    double min = std::numeric_limits<double>::max();
    double max = 0;
    const unsigned nnodes = check_allnodes ? element.getNumberOfNodes()
                                           : element.getNumberOfBaseNodes();
    for (unsigned i = 0; i < nnodes; i++)
    {
        for (unsigned j = i + 1; j < nnodes; j++)
        {
            const double dist(
                MathLib::sqrDist(*element.getNode(i), *element.getNode(j)));
            min = std::min(dist, min);
            max = std::max(dist, max);
        }
    }
    return {min, max};
}

std::pair<double, double> computeSqrEdgeLengthRange(Element const& element)
{
    double min = std::numeric_limits<double>::max();
    double max = 0;
    const unsigned nEdges(element.getNumberOfEdges());
    for (unsigned i = 0; i < nEdges; i++)
    {
        const double dist(MathLib::sqrDist(*element.getEdgeNode(i, 0),
                                           *element.getEdgeNode(i, 1)));
        min = std::min(dist, min);
        max = std::max(dist, max);
    }
    return {min, max};
}

bool isPointInElementXY(MathLib::Point3d const& p, Element const& e)
{
    for (std::size_t i(0); i < e.getNumberOfBaseNodes(); ++i)
    {
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

unsigned getNodeIDinElement(Element const& element, const MeshLib::Node* node)
{
    const unsigned nNodes(element.getNumberOfNodes());
    for (unsigned i(0); i < nNodes; i++)
    {
        if (node == element.getNode(i))
        {
            return i;
        }
    }
    return std::numeric_limits<unsigned>::max();
}

std::size_t getNodeIndex(Element const& element, unsigned const idx)
{
#ifndef NDEBUG
    if (idx >= element.getNumberOfNodes())
    {
        ERR("Error in MeshLib::getNodeIndex() - Index does not "
            "exist.");
        return std::numeric_limits<std::size_t>::max();
    }
#endif
    return element.getNode(idx)->getID();
}

}  // namespace MeshLib
