/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <vector>

#include "CoordinatesMappingTestData/TestHex8.h"
#include "CoordinatesMappingTestData/TestLine2.h"
#include "CoordinatesMappingTestData/TestLine3.h"
#include "CoordinatesMappingTestData/TestQuad4.h"
#include "CoordinatesMappingTestData/TestTri3.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalCoordinatesMapping.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "Tests/TestTools.h"

using namespace NumLib;
using namespace CoordinatesMappingTestData;

namespace
{
// Element types to be tested
using TestElementTypes =
    ::testing::Types<TestLine2, TestLine3, TestTri3, TestQuad4, TestHex8>;
}  // namespace

template <class T_TEST>
class NumLibFemNaturalCoordinatesMappingTest : public ::testing::Test,
                                               public T_TEST
{
public:
    using ElementType = typename T_TEST::ElementType;
    using ShapeFunctionType = typename T_TEST::ShapeFunctionType;
    static const unsigned dim = T_TEST::dim;
    static const unsigned e_nnodes = T_TEST::e_nnodes;
    static const unsigned global_dim = T_TEST::global_dim;
    // Matrix types
    using NodalVector = typename ::detail::EigenMatrixType<1, e_nnodes>::type;
    using DimNodalMatrix =
        typename ::detail::EigenMatrixType<dim, e_nnodes>::type;
    using DimMatrix = typename ::detail::EigenMatrixType<dim, dim>::type;
    using GlobalDimNodalMatrix =
        typename ::detail::EigenMatrixType<global_dim, e_nnodes>::type;
    // Shape data type
    using ShapeMatricesType = ShapeMatrices<NodalVector, DimNodalMatrix,
                                            DimMatrix, GlobalDimNodalMatrix>;
    // Natural coordinates mapping type
    using NaturalCoordsMappingType =
        NaturalCoordinatesMapping<ShapeFunctionType, ShapeMatricesType>;

public:
    NumLibFemNaturalCoordinatesMappingTest()
        : eps(2 * std::numeric_limits<double>::epsilon())
    {
        // create four elements used for testing
        naturalEle = this->createNaturalShape();
        irregularEle = this->createIrregularShape();
        clockwiseEle = this->createClockWise();
        zeroVolumeEle = this->createZeroVolume();

        // for destructor
        vec_eles.push_back(naturalEle);
        vec_eles.push_back(irregularEle);
        vec_eles.push_back(clockwiseEle);
        vec_eles.push_back(zeroVolumeEle);
        for (auto e : vec_eles)
        {
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                vec_nodes.push_back(e->getNode(i));
            }
        }
    }

    ~NumLibFemNaturalCoordinatesMappingTest() override
    {
        for (auto itr = vec_nodes.begin(); itr != vec_nodes.end(); ++itr)
        {
            delete *itr;
        }
        for (auto itr = vec_eles.begin(); itr != vec_eles.end(); ++itr)
        {
            delete *itr;
        }
    }

    const double eps;
    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const ElementType*> vec_eles;
    ElementType* naturalEle;
    ElementType* irregularEle;
    ElementType* clockwiseEle;
    ElementType* zeroVolumeEle;
};

TYPED_TEST_SUITE(NumLibFemNaturalCoordinatesMappingTest, TestElementTypes);

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_N)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    // only N
    NaturalCoordsMappingType::template computeShapeMatrices<ShapeMatrixType::N>(
        *this->naturalEle, this->r, shape, this->global_dim);
    ASSERT_FALSE(shape.N.isZero());
    ASSERT_TRUE(shape.dNdr.isZero());
    ASSERT_TRUE(shape.J.isZero());
    ASSERT_TRUE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_DNDR)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    // dNdr
    NaturalCoordsMappingType::template computeShapeMatrices<
        ShapeMatrixType::DNDR>(*this->naturalEle, this->r, shape,
                               this->global_dim);
    ASSERT_TRUE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_TRUE(shape.J.isZero());
    ASSERT_TRUE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_N_J)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    // N_J
    shape.setZero();
    NaturalCoordsMappingType::template computeShapeMatrices<
        ShapeMatrixType::N_J>(*this->naturalEle, this->r, shape,
                              this->global_dim);
    ASSERT_FALSE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest,
           CheckFieldSpecification_DNDR_J)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    // dNdr, J
    NaturalCoordsMappingType::template computeShapeMatrices<
        ShapeMatrixType::DNDR_J>(*this->naturalEle, this->r, shape,
                                 this->global_dim);
    ASSERT_TRUE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_DNDX)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    // DNDX
    shape.setZero();
    NaturalCoordsMappingType::template computeShapeMatrices<
        ShapeMatrixType::DNDX>(*this->naturalEle, this->r, shape, this->dim);
    ASSERT_TRUE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_FALSE(shape.invJ.isZero());
    ASSERT_FALSE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_ALL)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    // ALL
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices(*this->naturalEle, this->r,
                                                   shape, this->global_dim);
    ASSERT_FALSE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_FALSE(shape.invJ.isZero());
    ASSERT_FALSE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckNaturalShape)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;

    // identical to natural coordinates
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);
    NaturalCoordsMappingType::computeShapeMatrices(*this->naturalEle, this->r,
                                                   shape, this->global_dim);
    double exp_J[TestFixture::dim * TestFixture::dim] = {0.0};
    for (unsigned i = 0; i < this->dim; i++)
    {
        exp_J[i + this->dim * i] = 1.0;
    }

    ASSERT_ARRAY_NEAR(this->nat_exp_N, shape.N.data(), shape.N.size(),
                      this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdr.data(), shape.dNdr.size(),
                      this->eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), this->eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.invJ.data(), shape.invJ.size(), this->eps);
    ASSERT_NEAR(1.0, shape.detJ, this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdx.data(), shape.dNdx.size(),
                      this->eps);
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckIrregularShape)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;

    // irregular shape
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);
    NaturalCoordsMappingType::computeShapeMatrices(*this->irregularEle, this->r,
                                                   shape, this->global_dim);
    // std::cout <<  std::setprecision(16) << shape;

    ASSERT_ARRAY_NEAR(this->nat_exp_N, shape.N.data(), shape.N.size(),
                      this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdr.data(), shape.dNdr.size(),
                      this->eps);
    ASSERT_ARRAY_NEAR(this->ir_exp_J, shape.J.data(), shape.J.size(),
                      this->eps);
    ASSERT_NEAR(this->ir_exp_detJ, shape.detJ, this->eps);
    ASSERT_ARRAY_NEAR(this->ir_exp_invJ, shape.invJ.data(), shape.invJ.size(),
                      this->eps);
    ASSERT_ARRAY_NEAR(this->ir_exp_dNdx, shape.dNdx.data(), shape.dNdx.size(),
                      this->eps);
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckClockwise)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;

    // clockwise node ordering, which is invalid)
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    EXPECT_ANY_THROW(NaturalCoordsMappingType::computeShapeMatrices(
        *this->clockwiseEle, this->r, shape, this->global_dim));
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckZeroVolume)
{
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;
    using NaturalCoordsMappingType =
        typename TestFixture::NaturalCoordsMappingType;

    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);

    EXPECT_ANY_THROW(NaturalCoordsMappingType::computeShapeMatrices(
        *this->zeroVolumeEle, this->r, shape, this->global_dim));
}

TEST(NumLib, FemNaturalCoordinatesMappingLineY)
{
    using NodalVector = ::detail::EigenMatrixType<1, 2>::type;
    using DimNodalMatrix = ::detail::EigenMatrixType<1, 2>::type;
    using DimMatrix = ::detail::EigenMatrixType<1, 1>::type;
    using GlobalDimNodalMatrix = ::detail::EigenMatrixType<2, 2>::type;
    // Shape data type
    using ShapeMatricesType = ShapeMatrices<NodalVector, DimNodalMatrix,
                                            DimMatrix, GlobalDimNodalMatrix>;
    using MappingType =
        NaturalCoordinatesMapping<ShapeLine2, ShapeMatricesType>;
    double r[] = {0.5};
    auto line = TestLine2::createY();
    static const unsigned dim = 1;
    static const unsigned e_nnodes = 2;
    ShapeMatricesType shape(dim, 2, e_nnodes);
    MappingType::computeShapeMatrices(*line, r, shape, 2);

    double exp_J[dim * dim] = {0.0};
    for (unsigned i = 0; i < dim; i++)
    {
        exp_J[i + dim * i] = 1.0;
    }

    const double eps(std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(TestLine2::nat_exp_N, shape.N.data(), shape.N.size(),
                      eps);
    ASSERT_ARRAY_NEAR(TestLine2::nat_exp_dNdr, shape.dNdr.data(),
                      shape.dNdr.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.invJ.data(), shape.invJ.size(), eps);
    ASSERT_NEAR(1.0, shape.detJ, eps);
    double exp_dNdx[2 * e_nnodes] = {0, 0, -0.5, 0.5};
    ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), eps);

    for (auto n = 0u; n < line->getNumberOfNodes(); ++n)
    {
        delete line->getNode(n);
    }
    delete line;
}
