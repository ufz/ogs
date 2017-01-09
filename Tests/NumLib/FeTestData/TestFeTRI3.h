/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TESTFETRI3_H_
#define TESTFETRI3_H_

#include "MeshLib/Elements/Tri.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "MatrixTools.h"

namespace FeTestData
{

// Test case for TRI3
class TestFeTRI3
{
public:
    using ShapeFunction = NumLib::ShapeTri3;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeTRI3<ShapeMatrixPolicy_>;

    typedef MeshLib::Tri MeshElementType;
    static const unsigned dim = 2; //MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 3;
    static const unsigned n_sample_pt_order3 = 4;
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 0.0,  0.0,  0.0);
        nodes[1] = new MeshLib::Node( 1.0,  0.0,  0.0);
        nodes[2] = new MeshLib::Node( 0.0, 1.0,  0.0);
        return new MeshLib::Tri(nodes);
    }

    /// set an expected mass matrix for 1m x 1m
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.0; m(0,1) = 1./2; m(0,2) = 1./2;
        m(1,1) = 1.0; m(1,2) = 1./2;
        m(2,2) = 1.0;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= 1./6.*0.5;
    }

    /// set an expected laplace matrix for 1m x 1m
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 2.0; m(0,1) = -1.0; m(0,2) = -1.0;
        m(1,1) = 1.0; m(1,2) =  0.0;
        m(2,2) = 1.0;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= k*0.5;
    }
};

} // namespace

#endif

