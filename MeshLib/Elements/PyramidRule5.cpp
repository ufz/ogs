/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PyramidRule5.h"

#include <logog/include/logog.hpp>

#include "MathLib/GeometricBasics.h"

#include "MeshLib/Node.h"
#include "Quad.h"
#include "Tri.h"

namespace MeshLib {

const unsigned PyramidRule5::face_nodes[5][4] =
{
    {0, 1, 4, 99}, // Face 0
    {1, 2, 4, 99}, // Face 1
    {2, 3, 4, 99}, // Face 2
    {3, 0, 4, 99}, // Face 3
    {0, 3, 2,  1}  // Face 4
};

const unsigned PyramidRule5::edge_nodes[8][2] =
{
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {2, 3}, // Edge 2
    {0, 3}, // Edge 3
    {0, 4}, // Edge 4
    {1, 4}, // Edge 5
    {2, 4}, // Edge 6
    {3, 4}  // Edge 7
};

const unsigned PyramidRule5::n_face_nodes[5] = { 3, 3, 3, 3, 4 };

const Element* PyramidRule5::getFace(const Element* e, unsigned i)
{
    if (i<e->getNumberOfFaces())
    {
        unsigned nFaceNodes(PyramidRule5::n_face_nodes[i]);
        auto** nodes = new Node*[nFaceNodes];
        for (unsigned j=0; j<nFaceNodes; j++)
            nodes[j] = const_cast<Node*>(e->getNode(face_nodes[i][j]));

        if (i < 4)
            return new Tri(nodes);

        return new Quad(nodes);
    }
    ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
    return nullptr;
}

double PyramidRule5::computeVolume(Node const* const* _nodes)
{
    return MathLib::calcTetrahedronVolume(*_nodes[0], *_nodes[1], *_nodes[2], *_nodes[4])
         + MathLib::calcTetrahedronVolume(*_nodes[2], *_nodes[3], *_nodes[0], *_nodes[4]);
}

bool PyramidRule5::isPntInElement(Node const* const* nodes,
                                  MathLib::Point3d const& pnt,
                                  double eps)
{
    return (MathLib::isPointInTetrahedron(
                pnt, *nodes[0], *nodes[1], *nodes[2], *nodes[4], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[0], *nodes[2], *nodes[3], *nodes[4], eps));
}

unsigned PyramidRule5::identifyFace(Node const* const* _nodes, Node* nodes[3])
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

ElementErrorCode PyramidRule5::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();

    const auto* base(dynamic_cast<const MeshLib::Quad*>(e->getFace(4)));
    if (base)
    {
        error_code |= base->validate();
        error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    }
    else
        error_code.set(ElementErrorFlag::NodeOrder);
    delete base;

    return error_code;
}

} // end namespace MeshLib
