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

#include "NumLib/Fem/FemEnums.h"
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
    ASSERT_EQ(e_nnodes, shape.N.size());
    ASSERT_EQ(e_nnodes, shape.dNdr.cols());
    ASSERT_EQ(dim, shape.dNdr.rows());
    ASSERT_EQ(e_nnodes, shape.dNdx.cols());
    ASSERT_EQ(dim, shape.dNdx.rows());
    ASSERT_EQ(dim, shape.J.rows());
    ASSERT_EQ(dim, shape.J.cols());
    ASSERT_EQ(dim, shape.invJ.rows());
    ASSERT_EQ(dim, shape.invJ.cols());

    // check initialized with zero
    ASSERT_EQ(NodalVector::Zero(), shape.N);
    ASSERT_EQ(DimNodalMatrix::Zero(), shape.dNdr);
    ASSERT_EQ(DimNodalMatrix::Zero(), shape.dNdx);
    ASSERT_EQ(DimMatrix::Zero(), shape.J);
    ASSERT_EQ(DimMatrix::Zero(), shape.invJ);
    ASSERT_EQ(0.0, shape.detJ);

    // check setZero() for all
    setShapeDataToOnes(shape);
    shape.setZero();
    ASSERT_EQ(NodalVector::Zero(), shape.N);
    ASSERT_EQ(DimNodalMatrix::Zero(), shape.dNdr);
    ASSERT_EQ(DimNodalMatrix::Zero(), shape.dNdx);
    ASSERT_EQ(DimMatrix::Zero(), shape.J);
    ASSERT_EQ(DimMatrix::Zero(), shape.invJ);
    ASSERT_EQ(0.0, shape.detJ);

    // check setZero() only for N
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::N>();
    ASSERT_EQ(NodalVector::Zero(), shape.N);
    ASSERT_EQ(DimNodalMatrix::Ones(), shape.dNdr);
    ASSERT_EQ(DimNodalMatrix::Ones(), shape.dNdx);
    ASSERT_EQ(DimMatrix::Ones(), shape.J);
    ASSERT_EQ(DimMatrix::Ones(), shape.invJ);
    ASSERT_EQ(1.0, shape.detJ);

    // check setZero() for N, dNdr, and J
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::N_J>();
    ASSERT_EQ(NodalVector::Zero(), shape.N);
    ASSERT_EQ(DimNodalMatrix::Zero(), shape.dNdr);
    ASSERT_EQ(DimNodalMatrix::Ones(), shape.dNdx);
    ASSERT_EQ(DimMatrix::Zero(), shape.J);
    ASSERT_EQ(DimMatrix::Ones(), shape.invJ);
    ASSERT_EQ(0.0, shape.detJ);

    // check setZero() only for DNDX stuff
    setShapeDataToOnes(shape);
    shape.setZero<ShapeMatrixType::DNDX>();
    ASSERT_EQ(NodalVector::Ones(), shape.N);
    ASSERT_EQ(DimNodalMatrix::Zero(), shape.dNdr);
    ASSERT_EQ(DimNodalMatrix::Zero(), shape.dNdx);
    ASSERT_EQ(DimMatrix::Zero(), shape.J);
    ASSERT_EQ(DimMatrix::Zero(), shape.invJ);
    ASSERT_EQ(0.0, shape.detJ);
}
#endif


