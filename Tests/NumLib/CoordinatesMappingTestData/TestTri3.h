/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef COORDINATESMAPPINGTESTDATA_TESTTRI3_H_
#define COORDINATESMAPPINGTESTDATA_TESTTRI3_H_

#include "MeshLib/Elements/Tri.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"

namespace CoordinatesMappingTestData
{

class TestTri3
{
 public:
    // Element information
    typedef MeshLib::Tri ElementType;
    typedef NumLib::ShapeTri3 ShapeFunctionType;
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

    // element having shape identical to that in natural coordinates
    MeshLib::Tri* createNaturalShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 0.0,  0.0,  0.0);
        nodes[1] = new MeshLib::Node( 1.0,  0.0,  0.0);
        nodes[2] = new MeshLib::Node( 0.0, 1.0,  0.0);
        return new MeshLib::Tri(nodes);
    }

    // element having irregular or skew shape
    MeshLib::Tri* createIrregularShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 0.1, 0.1,  0.0);
        nodes[1] = new MeshLib::Node( 2.0, -0.6,  0.0);
        nodes[2] = new MeshLib::Node( 1.5,  0.4,  0.0);
        return new MeshLib::Tri(nodes);
    }

    // invalid case: clock wise node ordering
    MeshLib::Tri* createClockWise()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 0.0,  0.0,  0.0);
        nodes[2] = new MeshLib::Node( 1.0,  0.0,  0.0);
        nodes[1] = new MeshLib::Node( 0.0, 1.0,  0.0);
        return new MeshLib::Tri(nodes);
    }

    // invalid case: zero area
    MeshLib::Tri* createZeroVolume()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 0.0,  0.0,  0.0);
        nodes[1] = new MeshLib::Node( 1.0,  0.0,  0.0);
        nodes[2] = new MeshLib::Node( 1.0,  0.0,  0.0);
        return new MeshLib::Tri(nodes);
    }
};

const double TestTri3::r[dim] = {0.5, 0.5};
const double TestTri3::nat_exp_N[e_nnodes] = {0.0, 0.5, 0.5};
const double TestTri3::nat_exp_dNdr[e_nnodes*dim] = {-1, 1, 0, -1, 0, 1};
const double TestTri3::ir_exp_J[dim*dim] = {1.9, -0.7, 1.4, 0.3};
const double TestTri3::ir_exp_detJ = 1.55;
const double TestTri3::ir_exp_invJ[dim*dim]
    = {0.19354838709677424, 0.45161290322580649, -0.90322580645161299, 1.2258064516129032};
const double TestTri3::ir_exp_dNdx[dim*e_nnodes]
    = {-0.6451612903225807, 0.1935483870967742, 0.4516129032258065,
       -0.3225806451612903, -0.903225806451613, 1.225806451612903};
const double TestTri3::cl_exp_J[dim*dim] = {0.0, 1.0, 1.0, 0.0};
const double TestTri3::cl_exp_detJ = -1.;
const double TestTri3::ze_exp_J[dim*dim] = {1.0, 0.0, 1.0, 0.0};

}

#endif
