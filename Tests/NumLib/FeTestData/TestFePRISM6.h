/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Prism.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "MatrixTools.h"

namespace FeTestData
{

class TestFePRISM6
{
public:
    using ShapeFunction = NumLib::ShapePrism6;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FePRISM6<ShapeMatrixPolicy_>;

    typedef MeshLib::Prism MeshElementType;
    static const unsigned dim = 3;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 6;
    static const unsigned n_sample_pt_order3 = 6; //TODO no implementation yet
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        // cubic
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(0.5, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(0.5, 0.5, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 0.0, 0.5);
        nodes[4] = new MeshLib::Node(0.5, 0.0, 0.5);
        nodes[5] = new MeshLib::Node(0.5, 0.5, 0.5);
        return new MeshElementType(nodes);
    }

    /// set an expected mass matrix
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.; m(0,1) = 1./2.; m(0,2) = 1./2.; m(0,3) = 1./2.; m(0,4) = 1./4.; m(0,5) = 1./4.;
        m(1,1) = 1.; m(1,2) = 1./2.; m(1,3) = 1./4.; m(1,4) = 1./2.; m(1,5) = 1./4.;
        m(2,2) = 1.; m(2,3) = 1./4.; m(2,4) = 1./4.; m(2,5) = 1./2.;
        m(3,3) = 1.; m(3,4) = 1./2.; m(3,5) = 1./2.;
        m(4,4) = 1.; m(4,5) = 1./2.;
        m(5,5) = 1.;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= 1./288.;
    }

    /// set an expected laplace matrix
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.; m(0,1) = -1./2.; m(0,2) = 1./6.; m(0,3) = 0.0; m(0,4) = -1./2.; m(0,5) = -1./6.;
        m(1,1) = 5./3.; m(1,2) = -1./2.; m(1,3) = -1./2.; m(1,4) = 1./3.; m(1,5) = -1./2.;
        m(2,2) = 1.; m(2,3) = -1./6.; m(2,4) = -1./2.; m(2,5) = 0.0;
        m(3,3) = 1.; m(3,4) = -1./2.; m(3,5) = 1./6.;
        m(4,4) = 5./3.; m(4,5) = -1./2.;
        m(5,5) = 1.;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= k/8.;
    }
};

} // namespace
