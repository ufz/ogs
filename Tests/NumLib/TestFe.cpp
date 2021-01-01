/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <Eigen/Eigen>
#include <cmath>
#include <vector>

#include "FeTestData/TestFeHEX8.h"
#include "FeTestData/TestFeLINE2.h"
#include "FeTestData/TestFeLINE2Y.h"
#include "FeTestData/TestFeLINE3.h"
#include "FeTestData/TestFePRISM6.h"
#include "FeTestData/TestFePYRA5.h"
#include "FeTestData/TestFeQUAD4.h"
#include "FeTestData/TestFeTET4.h"
#include "FeTestData/TestFeTRI3.h"
#include "MeshLib/CoordinateSystem.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
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
template <class TestFeType_, template <typename, unsigned> class ShapeMatrixPolicy_>
struct TestCase
{
    using TestFeType = TestFeType_;
    static const unsigned GlobalDim = TestFeType::global_dim;
    using ShapeMatrixTypes = ShapeMatrixPolicy_<typename TestFeType::ShapeFunction, GlobalDim>;
    template <typename X>
    using ShapeMatrixPolicy = ShapeMatrixPolicy_<X, GlobalDim>;
    using ShapeFunction = typename TestFeType::ShapeFunction;
};

using TestTypes =
    ::testing::Types<TestCase<TestFeHEX8, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeLINE2, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeLINE2Y, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeLINE3, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFePRISM6, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFePYRA5, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeQUAD4, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeTET4, EigenDynamicShapeMatrixPolicy>,
                     TestCase<TestFeTRI3, EigenDynamicShapeMatrixPolicy>,

                     TestCase<TestFeHEX8, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeLINE2, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeLINE2Y, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeLINE3, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFePRISM6, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFePYRA5, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeQUAD4, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeTET4, EigenFixedShapeMatrixPolicy>,
                     TestCase<TestFeTRI3, EigenFixedShapeMatrixPolicy>>;
}  // namespace

template <class T>
class NumLibFemIsoTest : public ::testing::Test, public T::TestFeType
{
 public:
     using ShapeMatrixTypes = typename T::ShapeMatrixTypes;
     using ShapeFunction = typename T::ShapeFunction;
     using TestFeType = typename T::TestFeType;
     // Matrix types
     using NodalMatrix = typename ShapeMatrixTypes::NodalMatrixType;
     using NodalVector = typename ShapeMatrixTypes::NodalVectorType;
     using DimNodalMatrix = typename ShapeMatrixTypes::DimNodalMatrixType;
     using DimMatrix = typename ShapeMatrixTypes::DimMatrixType;
     using GlobalDimMatrixType = typename ShapeMatrixTypes::GlobalDimMatrixType;

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
     NumLibFemIsoTest()
         : D(dim, dim),
           expectedM(e_nnodes, e_nnodes),
           expectedK(e_nnodes, e_nnodes),
           integration_method(2)
     {
         // create a mesh element used for testing
         mesh_element = this->createMeshElement();

         // set a conductivity tensor
         setIdentityMatrix(dim, D);
         D *= conductivity;
         MeshLib::ElementCoordinatesMappingLocal ele_local_coord(
             *mesh_element,
             MeshLib::CoordinateSystem(*mesh_element).getDimension());
         auto R = ele_local_coord.getRotationMatrixToGlobal().topLeftCorner(
             TestFeType::dim, TestFeType::global_dim);
         globalD.noalias() = R.transpose() * (D * R);

         // set expected matrices
         this->setExpectedMassMatrix(expectedM);
         this->setExpectedLaplaceMatrix(conductivity, expectedK);

         // for destructor
         vec_eles.push_back(mesh_element);
         for (auto e : vec_eles)
         {
             for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
             {
                 vec_nodes.push_back(e->getNode(i));
             }
         }
    }

    ~NumLibFemIsoTest() override
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


    static const double conductivity;
    static const double eps;
    DimMatrix D;
    NodalMatrix expectedM;
    NodalMatrix expectedK;
    IntegrationMethod integration_method;
    GlobalDimMatrixType globalD;

    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const MeshElementType*> vec_eles;
    MeshElementType* mesh_element;

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
}; // NumLibFemIsoTest

template <class T>
const double NumLibFemIsoTest<T>::conductivity = 1e-11;

template <class T>
const double NumLibFemIsoTest<T>::eps = 2*std::numeric_limits<double>::epsilon();

template <class T>
const unsigned NumLibFemIsoTest<T>::dim;

template <class T>
const unsigned NumLibFemIsoTest<T>::e_nnodes;

template <class T>
const unsigned NumLibFemIsoTest<T>::n_sample_pt_order2;

template <class T>
const unsigned NumLibFemIsoTest<T>::n_sample_pt_order3;

TYPED_TEST_CASE(NumLibFemIsoTest, TestTypes);

TYPED_TEST(NumLibFemIsoTest, CheckMassMatrix)
{
    // Refer to typedefs in the fixture
    using NodalMatrix = typename TestFixture::NodalMatrix;

    auto const shape_matrices =
        NumLib::initShapeMatrices<typename TestFixture::ShapeFunction,
                                  typename TestFixture::ShapeMatrixTypes,
                                  TestFixture::global_dim,
                                  ShapeMatrixType::N_J>(
            *this->mesh_element, false /*is_axially_symmetric*/,
            this->integration_method);

    // evaluate a mass matrix M = int{ N^T D N }dA_e
    NodalMatrix M = NodalMatrix::Zero(this->e_nnodes, this->e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNumberOfPoints(); i++) {
        auto const& shape = shape_matrices[i];
        auto wp = this->integration_method.getWeightedPoint(i);
        M.noalias() += shape.N.transpose() * shape.N * shape.detJ * wp.getWeight();
    }
    //std::cout << "M=\n" << M;
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoTest, CheckLaplaceMatrix)
{
    // Refer to typedefs in the fixture
    using NodalMatrix = typename TestFixture::NodalMatrix;

    auto const shape_matrices =
        NumLib::initShapeMatrices<typename TestFixture::ShapeFunction,
                                  typename TestFixture::ShapeMatrixTypes,
                                  TestFixture::global_dim,
                                  ShapeMatrixType::DNDX>(
            *this->mesh_element, false /*is_axially_symmetric*/,
            this->integration_method);

    // evaluate a Laplace matrix K = int{ dNdx^T D dNdx }dA_e
    NodalMatrix K = NodalMatrix::Zero(this->e_nnodes, this->e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNumberOfPoints(); i++) {
        auto const& shape = shape_matrices[i];
        auto wp = this->integration_method.getWeightedPoint(i);
        K.noalias() += shape.dNdx.transpose() * this->globalD * shape.dNdx * shape.detJ * wp.getWeight();
    }
    //std::cout << "K=\n" << K << std::endl;
    ASSERT_ARRAY_NEAR(this->expectedK.data(), K.data(), K.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoTest, CheckMassLaplaceMatrices)
{
    // Refer to typedefs in the fixture
    using NodalMatrix = typename TestFixture::NodalMatrix;

    auto const shape_matrices =
        NumLib::initShapeMatrices<typename TestFixture::ShapeFunction,
                                  typename TestFixture::ShapeMatrixTypes,
                                  TestFixture::global_dim>(
            *this->mesh_element, false /*is_axially_symmetric*/,
            this->integration_method);

    // evaluate both mass and laplace matrices at once
    NodalMatrix M = NodalMatrix::Zero(this->e_nnodes, this->e_nnodes);
    NodalMatrix K = NodalMatrix::Zero(this->e_nnodes, this->e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNumberOfPoints(); i++) {
        auto const& shape = shape_matrices[i];
        auto wp = this->integration_method.getWeightedPoint(i);
        M.noalias() += shape.N.transpose() * shape.N * shape.detJ * wp.getWeight();
        K.noalias() += shape.dNdx.transpose() * (this->globalD * shape.dNdx) * shape.detJ * wp.getWeight();
    }

    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->expectedK.data(), K.data(), K.size(), this->eps);
}

#if 0
TYPED_TEST(NumLibFemIsoTest, CheckGaussLegendreIntegrationLevel)
{
    // Refer to typedefs in the fixture
    using NodalMatrix = typename TestFixture::NodalMatrix;

    auto shape_matrices =
        NumLib::initShapeMatrices<typename TestFixture::ShapeFunction,
                                  typename TestFixture::ShapeMatrixTypes,
                                  TestFixture::global_dim>(
            *this->mesh_element, false /*is_axially_symmetric*/,
            this->integration_method);

    // evaluate a mass matrix
    NodalMatrix M(this->e_nnodes, this->e_nnodes);
    M.setZero();
    ASSERT_EQ(TestFixture::n_sample_pt_order2, this->integration_method.getNumberOfPoints());
    for (std::size_t i=0; i < this->integration_method.getNumberOfPoints(); i++) {
        auto const& shape = shape_matrices[i];
        auto wp = this->integration_method.getWeightedPoint(i);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }
    //std::cout << "M=\n" << M << std::endl;
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);

    // Change gauss quadrature level to 3
    this->integration_method.setIntegrationOrder(3);
    M *= .0;
    ASSERT_EQ(TestFixture::n_sample_pt_order3, this->integration_method.getNumberOfPoints());

    shape_matrices =
        NumLib::initShapeMatrices<typename TestFixture::ShapeFunction,
                                  typename TestFixture::ShapeMatrixTypes,
                                  TestFixture::global_dim>(
            *this->mesh_element, false /*is_axially_symmetric*/,
            this->integration_method);

    for (std::size_t i=0; i < this->integration_method.getNumberOfPoints(); i++) {
        auto const& shape = shape_matrices[i];
        auto wp = this->integration_method.getWeightedPoint(i);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }
    //std::cout << "M=\n" << M << std::endl;
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}
#endif

