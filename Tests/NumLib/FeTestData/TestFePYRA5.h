/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TESTFEPYRA5_H_
#define TESTFEPYRA5_H_

#include "MeshLib/Elements/Pyramid.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "MatrixTools.h"

namespace FeTestData
{

class TestFePYRA5
{
public:
    using ShapeFunction = NumLib::ShapePyra5;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FePYRA5<ShapeMatrixPolicy_>;

    typedef MeshLib::Pyramid MeshElementType;
    static const unsigned dim = 3; //MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 5;
    static const unsigned n_sample_pt_order3 = 13;
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        // cubic
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.5);
        nodes[1] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(0.0, 0.5, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 0.5, 0.5);
        nodes[4] = new MeshLib::Node(0.25, 0.25, 0.25);
        return new MeshElementType(nodes);
    }

    /// set an expected mass matrix
    template <class T_MATRIX>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.390972222222e-03; m(0,1) = 6.923611111111e-04; m(0,2) = 3.493055555556e-04; m(0,3) = 6.923611111111e-04; m(0,4) = 125./16.*1e-4;
        m(1,1) = 1.390972222222e-03; m(1,2) = 6.923611111111e-04; m(1,3) = 3.493055555556e-04; m(1,4) = 125./16.*1e-4;
        m(2,2) = 1.390972222222e-03; m(2,3) = 6.923611111111e-04; m(2,4) = 125./16.*1e-4;
        m(3,3) = 1.390972222222e-03; m(3,4) = 125./16.*1e-4;
        m(4,4) = 1./480.;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
    }

    /// set an expected laplace matrix
    template <class T_MATRIX>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 7.654800000000e-2; m(0,1) = 6.785333333333e-3; m(0,2) = -6.785333333333e-3; m(0,3) = 6.785333333333e-3; m(0,4) = -8.333333333333e-2;
        m(1,1) = 7.654800000000e-2; m(1,2) = 6.785333333333e-3; m(1,3) = -6.785333333333e-3; m(1,4) = -8.333333333333e-2;
        m(2,2) = 7.654800000000e-2; m(2,3) = 6.785333333333e-3; m(2,4) = -8.333333333333e-2;
        m(3,3) = 7.654800000000e-2; m(3,4) = -8.333333333333e-2;
        m(4,4) = 1./3.;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= k;
    }
};

} // namespace

#endif

