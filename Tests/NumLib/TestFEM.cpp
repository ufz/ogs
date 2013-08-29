/**
 * \file TestFEM.cpp
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <cmath>
#include <Eigen>

#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "Tests/TestTools.h"

namespace
{

// set a expected laplace matrix for a square quad element (1m x 1m)
template <class T_MATRIX>
void setExpectedLaplaceMatrix(double k, T_MATRIX &_m)
{
    // set upper triangle entries
    _m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0;
    _m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
    _m(2,2) = 4.0; _m(2,3) = -1.0;
    _m(3,3) = 4.0;
    // copy the upper triangle to lower to make the matrix symmetric
    for (std::size_t i=0; i<4; i++)
        for (std::size_t j=0; j<i; j++)
            _m(i,j) = _m(j,i);
    //_m *= k / 6.0;
    for (std::size_t i=0; i<4; i++)
        for (std::size_t j=0; j<4; j++)
            _m(i,j) *= k/6.0;
}

} // namespace

TEST(NumLib, FEM_ElementIsoQuad4)
{
    // Local dense matrix type (row-majored), vector type
#define USE_FIXED_EIGEN
#ifdef USE_FIXED_EIGEN
    const static unsigned dim = 2;
    const static unsigned e_nnodes = 4;
    typedef Eigen::Matrix<double, e_nnodes, e_nnodes, Eigen::RowMajor> NodalMatrix;
    typedef Eigen::Matrix<double, dim, e_nnodes, Eigen::RowMajor> DimNodalMatrix;
    typedef Eigen::Matrix<double, dim, dim, Eigen::RowMajor> DimMatrix;
    typedef Eigen::Matrix<double, e_nnodes, 1> NodalVector;
#else
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> NodalMatrix;
    typedef NodalMatrix DimNodalMatrix;
    typedef NodalMatrix DimMatrix;
    typedef Eigen::VectorXd NodalVector;
#endif
    typedef Eigen::VectorXd LocalVector;

    // prepare a quad mesh element
    MeshLib::Node node1(0.0, 0.0, 0.0, 1);
    MeshLib::Node node2(1.0, 0.0, 0.0, 2);
    MeshLib::Node node3(1.0, 1.0, 0.0, 3);
    MeshLib::Node node4(0.0, 1.0, 0.0, 4);
    MeshLib::Node** nodes = new MeshLib::Node*[4];
    nodes[0] = &node1;
    nodes[1] = &node2;
    nodes[2] = &node3;
    nodes[3] = &node4;
    MeshLib::Quad e(nodes);

    // create a finite element object
#ifdef CXX_TEMPLATE_ALIASES_SUPPORTED
    typedef NumLib::QUAD4<NodalVector, DimNodalMatrix, DimMatrix> QUAD4Type;
#else
    typedef NumLib::QUAD4<NodalVector, DimNodalMatrix, DimMatrix>::type QUAD4Type;
#endif
    typedef typename QUAD4Type::ShapeDataType ShapeData;
    QUAD4Type fe(e, 2);

    // evaluate a local matrix K = int{ dN^T D dN }
    NodalMatrix K(e.getNNodes(), e.getNNodes());
    const double conductivity = 1e-11;
#ifdef USE_FIXED_EIGEN
    DimMatrix D = conductivity * DimMatrix::Identity();
#else
    DimMatrix D = conductivity * DimMatrix::Identity(e.getDimension(), e.getDimension());
#endif
    double x[3] = {};
    ShapeData shape(e.getDimension(), e.getNNodes());
    auto q = fe.getIntegrationMethod();
    for (size_t i=0; i<q.getNPoints(); i++) {
        double w = q.getPoint(i, x);
        fe.computeShapeFunctions(x, shape);
        K.noalias() += shape.dNdx.transpose() * D * shape.dNdx * shape.detJ * w;
    }
    NodalMatrix expectedK(e.getNNodes(), e.getNNodes());
    setExpectedLaplaceMatrix(conductivity, expectedK);
    ASSERT_DOUBLE_ARRAY_EQ(expectedK.data(), K.data(), K.rows()*K.cols(), 1e-12);

    // change the sampling level to 3
    q.setSamplingLevel(3);
    K.setZero(e.getNNodes(), e.getNNodes());
    for (size_t i=0; i<q.getNPoints(); i++) {
        double w = q.getPoint(i, x);
        fe.computeShapeFunctions(x, shape);
        K.noalias() += shape.dNdx.transpose() * D * shape.dNdx * shape.detJ * w;
    }
    ASSERT_DOUBLE_ARRAY_EQ(expectedK.data(), K.data(), K.rows()*K.cols(), 1e-12);

    // interpolate
    NodalVector vec_nodal_values1(e.getNNodes());
    vec_nodal_values1 << 0, 1, 1, 0;
    double x_node0[2] = {1,1};
    fe.computeShapeFunctions(x_node0, shape);
    ASSERT_NEAR(.0, fe.interpolate(shape, vec_nodal_values1), 1e-5);
    double x_node1[2] = {-1,1};
    fe.computeShapeFunctions(x_node1, shape);
    ASSERT_NEAR(1.0, fe.interpolate(shape, vec_nodal_values1), 1e-5);
    double x_node_mid[2] = {0,0};
    fe.computeShapeFunctions(x_node_mid, shape);
    ASSERT_NEAR(0.5, fe.interpolate(shape, vec_nodal_values1), 1e-5);

    // extrapolate
    q.setSamplingLevel(2);
    LocalVector vec_gp_values(q.getNPoints());
    for (std::size_t i=0; i<q.getNPoints(); i++) {
        q.getPoint(i, x);
        fe.computeShapeFunctions(x, shape);
        vec_gp_values(i) = fe.interpolate(shape, vec_nodal_values1);
    }
    NodalVector vec_nodal_values2;
    fe.extrapolate(vec_gp_values, vec_nodal_values2);
    ASSERT_DOUBLE_ARRAY_EQ(vec_nodal_values1, vec_nodal_values2, vec_nodal_values1.size(), 1e-12);
}

