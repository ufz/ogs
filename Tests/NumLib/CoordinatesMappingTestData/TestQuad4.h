/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"

namespace CoordinatesMappingTestData
{

class TestQuad4
{
 public:
    // Element information
    typedef MeshLib::Quad ElementType;
    typedef NumLib::ShapeQuad4 ShapeFunctionType;
    static const unsigned global_dim = ElementType::dimension;
    static const unsigned dim = 2; //ElementType::dimension;
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

    // element shape identical to that in natural coordinates
    MeshLib::Quad* createNaturalShape()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 1.0,  1.0,  0.0);
        nodes[1] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[2] = new MeshLib::Node(-1.0, -1.0,  0.0);
        nodes[3] = new MeshLib::Node( 1.0, -1.0,  0.0);
        return new MeshLib::Quad(nodes);
    }

    // element having irregular or skew shape
    MeshLib::Quad* createIrregularShape()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-0.5, -0.5,  0.0);
        nodes[1] = new MeshLib::Node( 0.6, -0.6,  0.0);
        nodes[2] = new MeshLib::Node( 0.5,  0.4,  0.0);
        nodes[3] = new MeshLib::Node(-0.3,  0.1,  0.0);
        return new MeshLib::Quad(nodes);
    }

    // invalid case: clock wise node ordering
    MeshLib::Quad* createClockWise()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 1.0,  1.0,  0.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[2] = new MeshLib::Node(-1.0, -1.0,  0.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0,  0.0);
        return new MeshLib::Quad(nodes);
    }

    // invalid case: zero area
    MeshLib::Quad* createZeroVolume()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 1.0,  1.0,  0.0);
        nodes[1] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[2] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[3] = new MeshLib::Node( 1.0,  1.0,  0.0);
        return new MeshLib::Quad(nodes);
    }
};

const double TestQuad4::r[dim] = {0.5, 0.5};
const double TestQuad4::nat_exp_N[e_nnodes] = {0.5625, 0.1875, 0.0625, 0.1875};
const double TestQuad4::nat_exp_dNdr[e_nnodes*dim]
    = {0.375, -0.375, -0.125, 0.125, 0.375, 0.125, -0.125, -0.375};
const double TestQuad4::ir_exp_J[dim*dim] = {-0.5125, 0.0, -0.0625, -0.35};
const double TestQuad4::ir_exp_detJ = 0.179375;
const double TestQuad4::ir_exp_invJ[dim*dim]
    = {-1.9512195121951219, 0.0, 0.3484320557491290, -2.8571428571428572};
const double TestQuad4::ir_exp_dNdx[dim*e_nnodes]
    = {-0.73170731707317072, 0.73170731707317072, 0.243902439024390, -0.24390243902439029,
       -0.940766550522648, -0.48780487804878048, 0.313588850174216, 1.1149825783972125};
const double TestQuad4::cl_exp_J[dim*dim] = {0.0, 1.0, 1.0, 0.0};
const double TestQuad4::cl_exp_detJ = -1.;
const double TestQuad4::ze_exp_J[dim*dim] = {1.0, 0.0, 0.0, 0.0};

}
