/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef COORDINATESMAPPINGTESTDATA_TESTHEX8_H_
#define COORDINATESMAPPINGTESTDATA_TESTHEX8_H_

#include "MeshLib/Elements/Hex.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"

namespace CoordinatesMappingTestData
{

class TestHex8
{
 public:
    // Element information
    typedef MeshLib::Hex ElementType;
    typedef NumLib::ShapeHex8 ShapeFunctionType;
    static const unsigned global_dim = ElementType::dimension;
    static const unsigned dim = 3; //ElementType::dimension;
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

    // element shape identical to that in natural coordinates (see ShapeHex8.h)
    MeshLib::Hex* createNaturalShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0, -1.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0, -1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0, -1.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0, -1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[5] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[7] = new MeshLib::Node(-1.0,  1.0,  1.0);
        return new MeshLib::Hex(nodes);
    }

    // element having irregular or skew shape
    MeshLib::Hex* createIrregularShape()
    {
        // two times longer in z direction than the natural
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0, -1.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0, -1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0, -1.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0, -1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  3.0);
        nodes[5] = new MeshLib::Node( 1.0, -1.0,  3.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  3.0);
        nodes[7] = new MeshLib::Node(-1.0,  1.0,  3.0);
        return new MeshLib::Hex(nodes);
    }

    // invalid case: clock wise node ordering
    MeshLib::Hex* createClockWise()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0, -1.0);
        nodes[3] = new MeshLib::Node( 1.0, -1.0, -1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0, -1.0);
        nodes[1] = new MeshLib::Node(-1.0,  1.0, -1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[7] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[5] = new MeshLib::Node(-1.0,  1.0,  1.0);
        return new MeshLib::Hex(nodes);
    }

    // invalid case: zero volume
    MeshLib::Hex* createZeroVolume()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0,  1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[5] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[7] = new MeshLib::Node(-1.0,  1.0,  1.0);
        return new MeshLib::Hex(nodes);
    }
};

const double TestHex8::r[dim] = {0.5, 0.5, 0.5};
const double TestHex8::nat_exp_N[e_nnodes]
    = {0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625};
const double TestHex8::nat_exp_dNdr[e_nnodes*dim]
    = {-0.03125, 0.03125, 0.09375, -0.09375, -0.09375, 0.09375, 0.28125, -0.28125,
       -0.03125, -0.09375, 0.09375, 0.03125, -0.09375, -0.28125,  0.28125,  0.09375,
       -0.03125, -0.09375, -0.28125, -0.09375,  0.03125,  0.09375,  0.28125,  0.09375};
const double TestHex8::ir_exp_J[dim*dim]= {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0};
const double TestHex8::ir_exp_detJ = 2.;
const double TestHex8::ir_exp_invJ[dim*dim] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1./2.0};
const double TestHex8::ir_exp_dNdx[dim*e_nnodes]
    = {-0.03125, 0.03125, 0.09375, -0.09375, -0.09375, 0.09375, 0.28125, -0.28125,
       -0.03125, -0.09375, 0.09375, 0.03125, -0.09375, -0.28125,  0.28125,  0.09375,
       -0.015625, -0.046875, -0.140625, -0.046875,  0.015625,  0.046875,  0.140625,  0.046875};
const double TestHex8::cl_exp_J[dim*dim] = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
const double TestHex8::cl_exp_detJ = -1.;
const double TestHex8::ze_exp_J[dim*dim] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

}

#endif //OGS_USE_EIGEN


