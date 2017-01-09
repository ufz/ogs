/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Line.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "MatrixTools.h"

namespace FeTestData
{

class TestFeLINE3
{
public:
    using ShapeFunction = NumLib::ShapeLine3;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeLINE3<ShapeMatrixPolicy_>;

    typedef MeshLib::Line3 MeshElementType;
    static const unsigned dim = MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 2;
    static const unsigned n_sample_pt_order3 = 3;
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(0.5, 0.0, 0.0);
        return new MeshElementType(nodes);
    }

    /// set an expected mass matrix
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1./9.; m(0,1) = -1./18.; m(0,2) = 1./9.;
        m(1,1) = 1./9.; m(1,2) = 1./9.;
        m(2,2) = 4./9.;
        // make symmetric
        copyUpperToLower(3, m);
    }

    /// set an expected laplace matrix
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 7./3.; m(0,1) = 1./3.; m(0,2) = -8./3.;
        m(1,1) = 7./3.; m(1,2) = -8./3.;
        m(2,2) = 16./3.;
        // make symmetric
        copyUpperToLower(3, m);
        m *= k;
    }
};

} // namespace
