/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <vector>
#include <cmath>
#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "FeTestData/TestFeLINE2.h"
#include "FeTestData/TestFeTRI3.h"
#include "FeTestData/TestFeQUAD4.h"
#include "FeTestData/TestFeHEX8.h"

#include "Tests/TestTools.h"

using namespace NumLib;
using namespace FeTestData;

namespace
{

// test cases
template <class TestFeType_, template <typename> class ShapeMatrixPolicy_>
struct TestCase
{
    typedef TestFeType_ TestFeType;
    typedef ShapeMatrixPolicy_<typename TestFeType::ShapeFunction> ShapeMatrixTypes;
};

typedef ::testing::Types<
#ifdef OGS_USE_EIGEN
        TestCase<TestFeLINE2, EigenDynamicShapeMatrixPolicy>,
        TestCase<TestFeTRI3, EigenDynamicShapeMatrixPolicy>,
        TestCase<TestFeQUAD4, EigenDynamicShapeMatrixPolicy>,
        TestCase<TestFeHEX8, EigenDynamicShapeMatrixPolicy>,
        TestCase<TestFeLINE2, EigenFixedShapeMatrixPolicy>,
        TestCase<TestFeTRI3, EigenFixedShapeMatrixPolicy>,
        TestCase<TestFeQUAD4, EigenFixedShapeMatrixPolicy>,
        TestCase<TestFeHEX8, EigenFixedShapeMatrixPolicy>
#endif
        > TestTypes;
}

template <class T>
class NumLibFemIsoTest : public ::testing::Test, public T::TestFeType
{
 public:
    typedef typename T::ShapeMatrixTypes ShapeMatrixTypes;
    typedef typename T::TestFeType TestFeType;
    // Matrix types
    typedef typename ShapeMatrixTypes::NodalMatrixType NodalMatrix;
    typedef typename ShapeMatrixTypes::NodalVectorType NodalVector;
    typedef typename ShapeMatrixTypes::DimNodalMatrixType DimNodalMatrix;
    typedef typename ShapeMatrixTypes::DimMatrixType DimMatrix;
    // Finite element type
    typedef typename TestFeType::template FeType<ShapeMatrixTypes>::type FeType;
    // Shape matrix data type
    typedef typename ShapeMatrixTypes::ShapeMatrices ShapeMatricesType;
    typedef typename TestFeType::MeshElementType MeshElementType;

    static const unsigned dim = TestFeType::dim;
    static const unsigned e_nnodes = TestFeType::e_nnodes;
    static const unsigned n_sample_pt_order2 = TestFeType::n_sample_pt_order2;
    static const unsigned n_sample_pt_order3 = TestFeType::n_sample_pt_order3;

    using IntegrationMethod =
        typename NumLib::GaussIntegrationPolicy<MeshElementType>::IntegrationMethod;


 public:
    NumLibFemIsoTest() :
        D(dim, dim),
        expectedM(e_nnodes,e_nnodes),
        expectedK(e_nnodes,e_nnodes),
        integration_method(2)
    {
        // create a mesh element used for testing
        mesh_element = this->createMeshElement();

        // set a conductivity tensor
        setIdentityMatrix(dim, D);
        D *= conductivity;

        // set expected matrices
        this->setExpectedMassMatrix(expectedM);
        this->setExpectedLaplaceMatrix(conductivity, expectedK);

        // for destructor
        vec_eles.push_back(mesh_element);
        for (auto e : vec_eles)
            for (unsigned i=0; i<e->getNNodes(true); i++)
                vec_nodes.push_back(e->getNode(i));
    }

    virtual ~NumLibFemIsoTest()
    {
        for (auto itr = vec_nodes.begin(); itr!=vec_nodes.end(); ++itr )
            delete *itr;
        for (auto itr = vec_eles.begin(); itr!=vec_eles.end(); ++itr )
            delete *itr;
    }


    static const double conductivity;
    static const double eps;
    DimMatrix D;
    NodalMatrix expectedM;
    NodalMatrix expectedK;
    IntegrationMethod integration_method;

    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const MeshElementType*> vec_eles;
    MeshElementType* mesh_element;

}; // NumLibFemIsoTest

template <class T>
const double NumLibFemIsoTest<T>::conductivity = 1e-11;

template <class T>
const double NumLibFemIsoTest<T>::eps = std::numeric_limits<double>::epsilon();

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
    typedef typename TestFixture::FeType FeType;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeType fe(*this->mesh_element);

    // evaluate a mass matrix M = int{ N^T D N }dA_e
    NodalMatrix M(this->e_nnodes, this->e_nnodes);
    M.setZero();
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.template computeShapeFunctions<ShapeMatrixType::N_J>(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }

    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoTest, CheckLaplaceMatrix)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeType FeType;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeType fe(*this->mesh_element);

    // evaluate a Laplace matrix K = int{ dNdx^T D dNdx }dA_e
    NodalMatrix K(this->e_nnodes, this->e_nnodes);
    K.setZero();
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.template computeShapeFunctions<ShapeMatrixType::DNDX>(wp.getCoords(), shape);
        K.noalias() += shape.dNdx.transpose() * this->D * shape.dNdx * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedK.data(), K.data(), K.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoTest, CheckMassLaplaceMatrices)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeType FeType;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeType fe(*this->mesh_element);

    // evaluate both mass and laplace matrices at once
    NodalMatrix M(this->e_nnodes, this->e_nnodes);
    M.setZero();
    NodalMatrix K(this->e_nnodes, this->e_nnodes);
    K.setZero();
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.computeShapeFunctions(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
        K.noalias() += shape.dNdx.transpose() * this->D * shape.dNdx * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->expectedK.data(), K.data(), K.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoTest, CheckGaussIntegrationLevel)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeType FeType;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object with gauss quadrature level 2
    FeType fe(*this->mesh_element);

    // evaluate a mass matrix
    NodalMatrix M(this->e_nnodes, this->e_nnodes);
    M.setZero();
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    ASSERT_EQ(TestFixture::n_sample_pt_order2, this->integration_method.getNPoints());
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.computeShapeFunctions(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);

    // Change gauss quadrature level to 3
    this->integration_method.setIntegrationOrder(3);
    M *= .0;
    ASSERT_EQ(TestFixture::n_sample_pt_order3, this->integration_method.getNPoints());
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.computeShapeFunctions(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}


