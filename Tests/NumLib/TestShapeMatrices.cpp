/**
 * \author Norihiro Watanabe
 * \date   2013-09-06
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#ifdef OGS_USE_EIGEN
#include <Eigen>
#endif

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

#include "../TestTools.h"

using namespace NumLib;


#ifdef OGS_USE_EIGEN
TEST(NumLib, FemShapeMatricesWithEigen)
{
    const static unsigned dim = 2;
    const static unsigned e_nnodes = 4;

    // Eigen matrix types
    typedef Eigen::Matrix<double, e_nnodes, 1> NodalVector;
    typedef Eigen::Matrix<double, dim, e_nnodes, Eigen::RowMajor> DimNodalMatrix;
    typedef Eigen::Matrix<double, dim, dim, Eigen::RowMajor> DimMatrix;

    // Shape data type
    typedef ShapeMatrices<NodalVector,DimNodalMatrix,DimMatrix> ShapeMatricesType;

    auto setShapeDataToOnes = [](ShapeMatricesType &shape)
            {
                shape.N.setOnes();
                shape.dNdr.setOnes();
                shape.dNdx.setOnes();
                shape.J.setOnes();
                shape.invJ.setOnes();
                shape.detJ = 1.0;
            };

    ShapeMatricesType shape(dim, e_nnodes);

    // check dimension of shape matrices
    EXPECT_EQ(e_nnodes, shape.N.size());
    EXPECT_EQ(e_nnodes, shape.dNdr.cols());
    EXPECT_EQ(dim, shape.dNdr.rows());
    EXPECT_EQ(e_nnodes, shape.dNdx.cols());
    EXPECT_EQ(dim, shape.dNdx.rows());
    EXPECT_EQ(dim, shape.J.rows());
    EXPECT_EQ(dim, shape.J.cols());
    EXPECT_EQ(dim, shape.invJ.rows());
    EXPECT_EQ(dim, shape.invJ.cols());

    // check initialized with zero
    EXPECT_TRUE(shape.N.isZero());
    EXPECT_TRUE(shape.dNdr.isZero());
    EXPECT_TRUE(shape.dNdx.isZero());
    EXPECT_TRUE(shape.J.isZero());
    EXPECT_TRUE(shape.invJ.isZero());
    EXPECT_EQ(0.0, shape.detJ);

    // check setZero() for all
    setShapeDataToOnes(shape);
    shape.setZero();
    EXPECT_TRUE(shape.N.isZero());
    EXPECT_TRUE(shape.dNdr.isZero());
    EXPECT_TRUE(shape.dNdx.isZero());
    EXPECT_TRUE(shape.J.isZero());
    EXPECT_TRUE(shape.invJ.isZero());
    EXPECT_EQ(0.0, shape.detJ);

    // check setZero() only for N
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::N>();
    EXPECT_TRUE(shape.N.isZero());
    EXPECT_TRUE(shape.dNdr.isOnes());
    EXPECT_TRUE(shape.dNdx.isOnes());
    EXPECT_TRUE(shape.J.isOnes());
    EXPECT_TRUE(shape.invJ.isOnes());
    EXPECT_EQ(1.0, shape.detJ);

    // check setZero() for dNdr
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::DNDR>();
    EXPECT_TRUE(shape.N.isOnes());
    EXPECT_TRUE(shape.dNdr.isZero());
    EXPECT_TRUE(shape.dNdx.isOnes());
    EXPECT_TRUE(shape.J.isOnes());
    EXPECT_TRUE(shape.invJ.isOnes());
    EXPECT_EQ(1.0, shape.detJ);

    // check setZero() for N, dNdr, and J
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::N_J>();
    EXPECT_TRUE(shape.N.isZero());
    EXPECT_TRUE(shape.dNdr.isZero());
    EXPECT_TRUE(shape.dNdx.isOnes());
    EXPECT_TRUE(shape.J.isZero());
    EXPECT_TRUE(shape.invJ.isOnes());
    EXPECT_EQ(0.0, shape.detJ);

    // check setZero() for dNdr, and J
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::DNDR_J>();
    EXPECT_TRUE(shape.N.isOnes());
    EXPECT_TRUE(shape.dNdr.isZero());
    EXPECT_TRUE(shape.dNdx.isOnes());
    EXPECT_TRUE(shape.J.isZero());
    EXPECT_TRUE(shape.invJ.isOnes());
    EXPECT_EQ(0.0, shape.detJ);

    // check setZero() only for DNDX stuff
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::DNDX>();
    EXPECT_TRUE(shape.N.isOnes());
    EXPECT_TRUE(shape.dNdr.isZero());
    EXPECT_TRUE(shape.dNdx.isZero());
    EXPECT_TRUE(shape.J.isZero());
    EXPECT_TRUE(shape.invJ.isZero());
    EXPECT_EQ(0.0, shape.detJ);
}
#endif


