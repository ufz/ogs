/**
 * \file
 * \brief Test the gradient of shape functions via the numerical volume
 *        integration.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on July 14, 2017, 10:08 AM
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <cmath>
#include <vector>

#include "FeTestData/TestFeHEX20.h"
#include "FeTestData/TestFeHEX8.h"
#include "FeTestData/TestFeLINE2.h"
#include "FeTestData/TestFeLINE2Y.h"
#include "FeTestData/TestFeLINE3.h"
#include "FeTestData/TestFePRISM15.h"
#include "FeTestData/TestFePRISM6.h"
#include "FeTestData/TestFePYRA13.h"
#include "FeTestData/TestFePYRA5.h"
#include "FeTestData/TestFeQUAD4.h"
#include "FeTestData/TestFeQUAD8.h"
#include "FeTestData/TestFeQUAD9.h"
#include "FeTestData/TestFeTET10.h"
#include "FeTestData/TestFeTET4.h"
#include "FeTestData/TestFeTRI3.h"
#include "FeTestData/TestFeTRI6.h"
#include "MeshLib/CoordinateSystem.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshToolsLib/ComputeElementVolumeNumerically.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "Tests/TestTools.h"

using namespace NumLib;
using namespace FeTestData;

namespace
{
// test cases
template <class TestFeType_, template <typename, int> class ShapeMatrixPolicy_>
struct TestCase
{
    using TestFeType = TestFeType_;
    static const int GlobalDim = TestFeType::global_dim;
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
}  // namespace

template <class T>
class GradShapeFunctionTest : public ::testing::Test, public T::TestFeType
{
public:
    using ShapeMatrixTypes = typename T::ShapeMatrixTypes;
    using TestFeType = typename T::TestFeType;

    // Finite element type
    template <typename X>
    using ShapeMatrixPolicy = typename T::template ShapeMatrixPolicy<X>;
    using MeshElementType = typename TestFeType::MeshElementType;

    static const int dim = TestFeType::dim;
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
        {
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                vec_nodes.push_back(e->getNode(i));
            }
        }
    }

    ~GradShapeFunctionTest() override
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

    IntegrationMethod integration_method;

    const double element_volume;
    MeshElementType* mesh_element;

    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const MeshElementType*> vec_eles;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};  // NumLibFemIsoTest

template <class T>
const int GradShapeFunctionTest<T>::dim;

template <class T>
const unsigned GradShapeFunctionTest<T>::e_nnodes;

template <class T>
const unsigned GradShapeFunctionTest<T>::n_sample_pt_order2;

template <class T>
const unsigned GradShapeFunctionTest<T>::n_sample_pt_order3;

TYPED_TEST_SUITE(GradShapeFunctionTest, TestTypes);

TYPED_TEST(GradShapeFunctionTest,
           CheckGradShapeFunctionByComputingElementVolume)
{
    auto const shape_matrices =
        NumLib::initShapeMatrices<typename TestFixture::ShapeFunction,
                                  typename TestFixture::ShapeMatrixTypes,
                                  TestFixture::global_dim,
                                  ShapeMatrixType::N_J>(
            *this->mesh_element, false /*is_axially_symmetric*/,
            this->integration_method);

    // Compute element volume numerically as V_e = int {1}dA_e
    double computed_element_volume = 0.;
    for (std::size_t i = 0; i < this->integration_method.getNumberOfPoints();
         i++)
    {
        auto const& shape = shape_matrices[i];
        auto wp = this->integration_method.getWeightedPoint(i);
        computed_element_volume += shape.detJ * wp.getWeight();
    }

    ASSERT_NEAR(this->element_volume, computed_element_volume, 1.e-11);

    ASSERT_NEAR(
        MeshToolsLib::computeElementVolumeNumerically(*this->mesh_element),
        computed_element_volume, 1.e-11);
}
