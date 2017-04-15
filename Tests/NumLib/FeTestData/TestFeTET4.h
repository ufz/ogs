/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Tet.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "MatrixTools.h"

namespace FeTestData
{

class TestFeTET4
{
public:
    using ShapeFunction = NumLib::ShapeTet4;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeTET4<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Tet;
    static const unsigned dim = 3;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 5;
    static const unsigned n_sample_pt_order3 = 15;
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(0.5, 0.5, 0.5);
        nodes[2] = new MeshLib::Node(0.5, 0.0, 0.5);
        nodes[3] = new MeshLib::Node(1.0, 0.5, 0.5);
        return new MeshElementType(nodes);
    }

    /// set an expected mass matrix
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.0; m(0,1) = 1./2; m(0,2) = 1./2; m(0,3) = 1./2;
        m(1,1) = 1.0; m(1,2) = 1./2; m(1,3) = 1./2;
        m(2,2) = 1.0; m(2,3) = 1./2;
        m(3,3) = 1.0;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= 1./480.;
    }

    /// set an expected laplace matrix
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1./3.; m(0,1) = -1./3.; m(0,2) = -1./3.; m(0,3) = 1./3.;
        m(1,1) = 1.; m(1,2) = 0; m(1,3) = -2./3.;
        m(2,2) = 2./3.; m(2,3) = -1./3.;
        m(3,3) = 2./3.;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= k*0.25;
    }
};

} // namespace
