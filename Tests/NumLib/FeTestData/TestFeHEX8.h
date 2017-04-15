/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Hex.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "MatrixTools.h"

namespace FeTestData
{

// Test case for HEX8
class TestFeHEX8
{
public:
    using ShapeFunction = NumLib::ShapeHex8;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeHEX8<ShapeMatrixPolicy_>;

    typedef MeshLib::Hex MeshElementType;
    static const unsigned dim = 3; //MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 2*2*2;
    static const unsigned n_sample_pt_order3 = 3*3*3;
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshLib::Hex* createMeshElement()
    {
        // cubic
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0);
        nodes[4] = new MeshLib::Node(0.0, 0.0, 1.0);
        nodes[5] = new MeshLib::Node(1.0, 0.0, 1.0);
        nodes[6] = new MeshLib::Node(1.0, 1.0, 1.0);
        nodes[7] = new MeshLib::Node(0.0, 1.0, 1.0);
        return new MeshLib::Hex(nodes);
    }

    /// set an expected mass matrix
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.0; m(0,1) = 1./2; m(0,2) = 1./4; m(0,3) = 1./2; m(0,4) = 1./2; m(0,5) = 1./4; m(0,6) = 1./8; m(0,7) = 1./4;
        m(1,1) = 1.0; m(1,2) = 1./2; m(1,3) = 1./4; m(1,4) = 1./4; m(1,5) = 1./2; m(1,6) = 1./4; m(1,7) = 1./8;
        m(2,2) = 1.0; m(2,3) = 1./2; m(2,4) = 1./8; m(2,5) = 1./4; m(2,6) = 1./2; m(2,7) = 1./4;
        m(3,3) = 1.0; m(3,4) = 1./4; m(3,5) = 1./8; m(3,6) = 1./4; m(3,7) = 1./2;
        m(4,4) = 1.0; m(4,5) = 1./2; m(4,6) = 1./4; m(4,7) = 1./2;
        m(5,5) = 1.0; m(5,6) = 1./2; m(5,7) = 1./4;
        m(6,6) = 1.0; m(6,7) = 1./2;
        m(7,7) = 1.0;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= 1./27.;
    }

    /// set an expected laplace matrix
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 2./3; m(0,1) = 0; m(0,2) = -1./6; m(0,3) = 0; m(0,4) = 0; m(0,5) = -1./6; m(0,6) = -1./6; m(0,7) = -1./6;
        m(1,1) = 2./3; m(1,2) = 0; m(1,3) = -1./6; m(1,4) = -1./6; m(1,5) = 0; m(1,6) = -1./6; m(1,7) = -1./6;
        m(2,2) = 2./3; m(2,3) = 0; m(2,4) = -1./6; m(2,5) = -1./6; m(2,6) = 0; m(2,7) = -1./6;
        m(3,3) = 2./3; m(3,4) = -1./6; m(3,5) = -1./6; m(3,6) = -1./6; m(3,7) = 0;
        m(4,4) = 2./3; m(4,5) = 0; m(4,6) = -1./6; m(4,7) = 0;
        m(5,5) = 2./3; m(5,6) = 0; m(5,7) = -1./6;
        m(6,6) = 2./3; m(6,7) = 0;
        m(7,7) = 2./3;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= k/2.;
    }
};

} // namespace
