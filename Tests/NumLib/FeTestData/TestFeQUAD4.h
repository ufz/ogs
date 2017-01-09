/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "MatrixTools.h"

namespace FeTestData
{

// Test case for QUAD4
class TestFeQUAD4
{
public:
    using ShapeFunction = NumLib::ShapeQuad4;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeQUAD4<ShapeMatrixPolicy_>;

    typedef MeshLib::Quad MeshElementType;
    static const unsigned dim = 2; //MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 2*2;
    static const unsigned n_sample_pt_order3 = 3*3;
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshLib::Quad* createMeshElement()
    {
        // square
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0);
        return new MeshLib::Quad(nodes);
    }

    /// set an expected mass matrix for 1m x 1m
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.0; m(0,1) = 1./2; m(0,2) = 1./4; m(0,3) = 1./2;
        m(1,1) = 1.0; m(1,2) = 1./2; m(1,3) = 1./4;
        m(2,2) = 1.0; m(2,3) = 1./2;
        m(3,3) = 1.0;
        // make symmetric
        copyUpperToLower(4, m);
        m *= 1./9.;
    }

    /// set an expected laplace matrix for 1m x 1m
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 4.0; m(0,1) = -1.0; m(0,2) = -2.0; m(0,3) = -1.0;
        m(1,1) = 4.0; m(1,2) = -1.0; m(1,3) = -2.0;
        m(2,2) = 4.0; m(2,3) = -1.0;
        m(3,3) = 4.0;
        // make symmetric
        copyUpperToLower(4, m);
        m *= k/6.;
    }
};

} // namespace
