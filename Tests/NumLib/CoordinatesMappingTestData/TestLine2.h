/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include "MeshLib/Elements/Line.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"

namespace CoordinatesMappingTestData
{

class TestLine2
{
public:
    // Element information
    typedef MeshLib::Line ElementType;
    typedef NumLib::ShapeLine2 ShapeFunctionType;
    static const unsigned global_dim = ElementType::dimension;
    static const unsigned dim = ElementType::dimension;
    static const unsigned e_nnodes = ElementType::n_all_nodes;
    // Coordinates where shape functions are evaluated
    static const double r[dim];
    // Expected results for natural shape
    static const double nat_exp_N[e_nnodes];
    static const double nat_exp_dNdr[e_nnodes*dim];
    // Expected results for irregular shape
    static const double ir_exp_J[dim*dim];
    static const double ir_exp_invJ[dim*dim];
    static const double ir_exp_detJ;
    static const double ir_exp_dNdx[e_nnodes*dim];
    // Expected results for clock-wise node ordering
    static const double cl_exp_J[dim*dim];
    // Expected results for zero volume
    static const double cl_exp_detJ;
    static const double ze_exp_J[dim*dim];

    // element shape identical to that in natural coordinates (see ShapeLine2.h)
    static MeshLib::Line* createNaturalShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node( 1.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // element having irregular or skew shape
    MeshLib::Line* createIrregularShape()
    {
        // two times longer than the natural
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-2.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node( 2.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // invalid case: clock wise node ordering
    MeshLib::Line* createClockWise()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[1] = new MeshLib::Node(-1.0, 0.0, 0.0);
        nodes[0] = new MeshLib::Node( 1.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // invalid case: zero volume
    MeshLib::Line* createZeroVolume()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(0.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // 1.5d line
    static MeshLib::Line* createY()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, -1.0, 0.0);
        nodes[1] = new MeshLib::Node(0.0,  1.0, 0.0);
        return new MeshLib::Line(nodes);
    }
};

const double TestLine2::r[dim] = {0.5};
const double TestLine2::nat_exp_N[e_nnodes] = {0.25, 0.75};
const double TestLine2::nat_exp_dNdr[e_nnodes*dim] = {-0.5, 0.5};
const double TestLine2::ir_exp_J[dim*dim] = {2.0};
const double TestLine2::ir_exp_invJ[dim*dim] = {0.5};
const double TestLine2::ir_exp_detJ = 2.;
const double TestLine2::ir_exp_dNdx[dim*e_nnodes] = {-0.25, 0.25};
const double TestLine2::cl_exp_J[dim*dim] = {-1.};
const double TestLine2::cl_exp_detJ = -1;
const double TestLine2::ze_exp_J[dim*dim] = {0.0};

}
