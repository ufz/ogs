/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PrismRule6.h"

#include <logog/include/logog.hpp>

#include "MathLib/GeometricBasics.h"

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib {

const unsigned PrismRule6::face_nodes[5][4] =
{
    {0, 2, 1, 99}, // Face 0
    {0, 1, 4,  3}, // Face 1
    {1, 2, 5,  4}, // Face 2
    {2, 0, 3,  5}, // Face 3
    {3, 4, 5, 99}  // Face 4
};

const unsigned PrismRule6::edge_nodes[9][2] =
{
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {0, 2}, // Edge 2
    {0, 3}, // Edge 3
    {1, 4}, // Edge 4
    {2, 5}, // Edge 5
    {3, 4}, // Edge 6
    {4, 5}, // Edge 7
    {3, 5}  // Edge 8
};

const unsigned PrismRule6::n_face_nodes[5] = { 3, 4, 4, 4, 3 };

const Element* PrismRule6::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        unsigned nFaceNodes(PrismRule6::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j=0; j<nFaceNodes; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));

        if (i==0 || i==4)
            return new Tri(nodes);
        else
            return new Quad(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

double PrismRule6::computeVolume(Node const* const* _nodes)
{
    return MathLib::calcTetrahedronVolume(*_nodes[0], *_nodes[1], *_nodes[2], *_nodes[3])
         + MathLib::calcTetrahedronVolume(*_nodes[1], *_nodes[4], *_nodes[2], *_nodes[3])
         + MathLib::calcTetrahedronVolume(*_nodes[2], *_nodes[4], *_nodes[5], *_nodes[3]);
}

bool PrismRule6::isPntInElement(Node const* const* nodes,
                                MathLib::Point3d const& pnt,
                                double eps)
{
    return (MathLib::isPointInTetrahedron(
                pnt, *nodes[0], *nodes[1], *nodes[2], *nodes[3], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[1], *nodes[4], *nodes[2], *nodes[3], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[2], *nodes[4], *nodes[5], *nodes[3], eps));
}

unsigned PrismRule6::identifyFace(Node const* const* _nodes, Node* nodes[3])
{
    for (unsigned i=0; i<5; i++)
    {
        unsigned flag(0);
        for (unsigned j=0; j<4; j++)
            for (unsigned k=0; k<3; k++)
                if (face_nodes[i][j] != 99 && _nodes[face_nodes[i][j]] == nodes[k])
                    flag++;
        if (flag==3)
            return i;
    }
    return std::numeric_limits<unsigned>::max();
}

ElementErrorCode PrismRule6::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();

    for (unsigned i=1; i<4; ++i)
    {
        const auto* quad(dynamic_cast<const MeshLib::Quad*>(e->getFace(i)));
        if (quad)
            error_code |= quad->validate();
        else
            error_code.set(ElementErrorFlag::NodeOrder);
        delete quad;
    }
    error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    return error_code;
}

} // end namespace MeshLib
