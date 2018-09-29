/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

class TestFeLINE2Y
{
public:
    using ShapeFunction = NumLib::ShapeLine2;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeLINE2<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Line;
    static const unsigned dim = MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 2;
    static const unsigned n_sample_pt_order3 = 3;
    static const unsigned global_dim = 2;

    /// create a mesh element
    MeshLib::Line* createMeshElement()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(0.0, 1.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    double getVolume() const { return 1.0; }

    /// set an expected mass matrix
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1./3.; m(0,1) = 1./6.;
        m(1,1) = 1./3.;
        // make symmetric
        copyUpperToLower(2, m);
    }

    /// set an expected laplace matrix
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.0; m(0,1) = -1.0;
        m(1,1) = 1.0;
        // make symmetric
        copyUpperToLower(2, m);
        m *= k;
    }
};

} // namespace
