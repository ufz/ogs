/**
 * \brief Test the gradient of shape functions via the numerical volume
 *        integration.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestGradShapeFunction.cpp
 *  Created on July 14, 2017, 10:08 AM
 */

#include <gtest/gtest.h>

#include <vector>
#include <cmath>

#include <Eigen/Eigen>

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/CoordinateSystem.h"

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "FeTestData/TestFeLINE2.h"
#include "FeTestData/TestFeLINE2Y.h"
#include "FeTestData/TestFeLINE3.h"
#include "FeTestData/TestFeTRI3.h"
#include "FeTestData/TestFeTRI6.h"
#include "FeTestData/TestFeQUAD4.h"
#include "FeTestData/TestFeQUAD8.h"
#include "FeTestData/TestFeQUAD9.h"
#include "FeTestData/TestFeHEX8.h"
#include "FeTestData/TestFeHEX20.h"
#include "FeTestData/TestFeTET4.h"
#include "FeTestData/TestFeTET10.h"
#include "FeTestData/TestFePRISM6.h"
#include "FeTestData/TestFePRISM15.h"
#include "FeTestData/TestFePYRA5.h"
#include "FeTestData/TestFePYRA13.h"

#include "Tests/TestTools.h"

using namespace NumLib;
using namespace FeTestData;

namespace
{
// test cases
template <class TestFeType_,
          template <typename, unsigned> class ShapeMatrixPolicy_>
struct TestCase
{
    using TestFeType = TestFeType_;
    static const unsigned GlobalDim = TestFeType::global_dim;
    using ShapeMatrixTypes =
        ShapeMatrixPolicy_<typename TestFeType::ShapeFunction, GlobalDim>;
    template <typename X>
    using ShapeMatrixPolicy = ShapeMatrixPolicy_<X, GlobalDim>;
};

using TestTypes =
    ::testing::Types<TestCase<TestFeHEX8, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeHEX20, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeLINE2, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeLINE2Y, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeLINE3, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFePRISM6, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFePRISM15, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFePYRA5, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFePYRA13, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeQUAD4, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeQUAD8, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeQUAD9, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeTET4, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeTET10, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeTRI3, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeTRI6, EigenDynamicShapeMatrixPolicy>,

                     TestCase<TestFeHEX8, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeHEX20, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeLINE2, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeLINE2Y, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeLINE3, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFePRISM6, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFePRISM15, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFePYRA5, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFePYRA13, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeQUAD4, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeQUAD8, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeQUAD9, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeTET4, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeTET10, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeTRI3, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeTRI6, EigenFixedShapeMatrixPolicy>>;
}

template <class T>
class GradShapeFunctionTest : public ::testing::Test, public T::TestFeType
{
public:
    using ShapeMatrixTypes = typename T::ShapeMatrixTypes;
    using TestFeType = typename T::TestFeType;

    // Finite element type
    template <typename X>
    using ShapeMatrixPolicy = typename T::template ShapeMatrixPolicy<X>;
    using FeType =
        typename TestFeType::template FeType<ShapeMatrixPolicy>::type;

    // Shape matrix data type
    using ShapeMatricesType = typename ShapeMatrixTypes::ShapeMatrices;
    using MeshElementType = typename TestFeType::MeshElementType;

    static const unsigned dim = TestFeType::dim;
    static const unsigned e_nnodes = TestFeType::e_nnodes;
    static const unsigned n_sample_pt_order2 = TestFeType::n_sample_pt_order2;
    static const unsigned n_sample_pt_order3 = TestFeType::n_sample_pt_order3;

    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        MeshElementType>::IntegrationMethod;

public:
    GradShapeFunctionTest()
        : integration_method(3),
          element_volume(this->getVolume()),
          mesh_element(this->createMeshElement())
    {
        // only for destructor because class element has nodes in pointer type.
        vec_eles.push_back(mesh_element);
        for (auto e : vec_eles)
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
                vec_nodes.push_back(e->getNode(i));
    }

    ~GradShapeFunctionTest() override
    {
        for (auto itr = vec_nodes.begin(); itr != vec_nodes.end(); ++itr)
            delete *itr;
        for (auto itr = vec_eles.begin(); itr != vec_eles.end(); ++itr)
            delete *itr;
    }

    IntegrationMethod integration_method;

    const double element_volume;
    MeshElementType* mesh_element;

    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const MeshElementType*> vec_eles;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};  // NumLibFemIsoTest

template <class T>
const unsigned GradShapeFunctionTest<T>::dim;

template <class T>
const unsigned GradShapeFunctionTest<T>::e_nnodes;

template <class T>
const unsigned GradShapeFunctionTest<T>::n_sample_pt_order2;

template <class T>
const unsigned GradShapeFunctionTest<T>::n_sample_pt_order3;

TYPED_TEST_CASE(GradShapeFunctionTest, TestTypes);

TYPED_TEST(GradShapeFunctionTest,
           CheckGradShapeFunctionByComputingElementVolume)
{
    // Refer to typedefs in the fixture
    using FeType = typename TestFixture::FeType;
    using ShapeMatricesType = typename TestFixture::ShapeMatricesType;

    // create a finite element object
    FeType fe(*this->mesh_element);

    // Compute element volume numerically as V_e = int {1}dA_e
    double computed_element_volume = 0.;
    ShapeMatricesType shape(this->dim, this->global_dim, this->e_nnodes);
    for (std::size_t i = 0; i < this->integration_method.getNumberOfPoints();
         i++)
    {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.template computeShapeFunctions<ShapeMatrixType::N_J>(
            wp.getCoords(), shape, this->global_dim, false);
        computed_element_volume += shape.detJ * wp.getWeight();
    }

    ASSERT_NEAR(this->element_volume, computed_element_volume, 1.e-11);
}
