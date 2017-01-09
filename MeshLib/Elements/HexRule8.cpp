/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HexRule8.h"

#include <array>

#include <logog/include/logog.hpp>

#include "MathLib/GeometricBasics.h"

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Line.h"

namespace MeshLib {

const unsigned HexRule8::face_nodes[6][4] =
{
    {0, 3, 2, 1}, // Face 0
    {0, 1, 5, 4}, // Face 1
    {1, 2, 6, 5}, // Face 2
    {2, 3, 7, 6}, // Face 3
    {3, 0, 4, 7}, // Face 4
    {4, 5, 6, 7}  // Face 5
};

const unsigned HexRule8::edge_nodes[12][2] =
{
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {4, 5}, // Edge 4
    {5, 6}, // Edge 5
    {6, 7}, // Edge 6
    {4, 7}, // Edge 7
    {0, 4}, // Edge 8
    {1, 5}, // Edge 9
    {2, 6}, // Edge 10
    {3, 7}  // Edge 11
};

const Element* HexRule8::getFace(const Element* e, unsigned i)
{
    if (i < n_faces)
    {
        std::array<Node*, 4> nodes;
        for (unsigned j=0; j<4; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));
        return new Quad(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

double HexRule8::computeVolume(Node const* const* _nodes)
{
    return MathLib::calcTetrahedronVolume(*_nodes[4], *_nodes[7], *_nodes[5], *_nodes[0])
         + MathLib::calcTetrahedronVolume(*_nodes[5], *_nodes[3], *_nodes[1], *_nodes[0])
         + MathLib::calcTetrahedronVolume(*_nodes[5], *_nodes[7], *_nodes[3], *_nodes[0])
         + MathLib::calcTetrahedronVolume(*_nodes[5], *_nodes[7], *_nodes[6], *_nodes[2])
         + MathLib::calcTetrahedronVolume(*_nodes[1], *_nodes[3], *_nodes[5], *_nodes[2])
         + MathLib::calcTetrahedronVolume(*_nodes[3], *_nodes[7], *_nodes[5], *_nodes[2]);
}

bool HexRule8::isPntInElement(Node const* const* nodes,
                              MathLib::Point3d const& pnt,
                              double eps)
{
    return (MathLib::isPointInTetrahedron(
                pnt, *nodes[4], *nodes[7], *nodes[5], *nodes[0], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[5], *nodes[3], *nodes[1], *nodes[0], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[5], *nodes[7], *nodes[3], *nodes[0], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[5], *nodes[7], *nodes[6], *nodes[2], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[1], *nodes[3], *nodes[5], *nodes[2], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[3], *nodes[7], *nodes[5], *nodes[2], eps));
}

unsigned HexRule8::identifyFace(Node const* const* _nodes, Node* nodes[3])
{
    for (unsigned i=0; i<6; i++)
    {
        unsigned flag(0);
        for (unsigned j=0; j<4; j++)
            for (unsigned k=0; k<3; k++)
                if (_nodes[face_nodes[i][j]] == nodes[k])
                    flag++;
        if (flag==3)
            return i;
    }
    return std::numeric_limits<unsigned>::max();
}

ElementErrorCode HexRule8::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();

    for (unsigned i=0; i<6; ++i)
    {
        if (error_code.all())
            break;

        const MeshLib::Element* quad (e->getFace(i));
        error_code |= quad->validate();
        delete quad;
    }
    error_code[ElementErrorFlag::NodeOrder]  = !e->testElementNodeOrder();
    return error_code;
}

} // end namespace MeshLib
