/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MeshLib/Elements/Line.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "Tests/TestTools.h"

using namespace NumLib;

namespace
{

// copy matrix entries in upper triangle to lower triangle
template <class T_MATRIX, typename ID_TYPE=signed>
void copyUpperToLower(const ID_TYPE dim, T_MATRIX &m)
{
    for (ID_TYPE i=0; i<dim; i++)
        for (ID_TYPE j=0; j<i; j++)
            m(i,j) = m(j,i);
}

// set an identity matrix
template <class T_MATRIX, typename ID_TYPE=signed>
void setIdentityMatrix(unsigned dim, T_MATRIX &m)
{
    for (unsigned i=0; i<dim; i++)
        for (unsigned j=0; j<dim; j++)
            m(i,j) = 0.0;
    for (unsigned i=0; i<dim; i++)
        m(i,i) = 1.0;
}

// Test case for LINE2
class TestFeLINE2
{
public:
    template <class T_MATRIX_TYPES>
    using FeType = NumLib::FeLINE2<typename T_MATRIX_TYPES::NodalVectorType, typename T_MATRIX_TYPES::DimNodalMatrixType, typename T_MATRIX_TYPES::DimMatrixType>;
    typedef MeshLib::Line MeshElementType;
    static const unsigned dim = MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;

    MeshLib::Line* createMeshElement()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // set an expected mass matrix
    template <class T_MATRIX, typename ID_TYPE=signed>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1./3.; m(0,1) = 1./6.;
        m(1,1) = 1./3.;
        // make symmetric
        copyUpperToLower(2, m);
    }

    // set an expected laplace matrix
    template <class T_MATRIX, typename ID_TYPE=signed>
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

// Test case for QUAD4
class TestFeQUAD4
{
public:
    template <class T_MATRIX_TYPES>
    using FeType = NumLib::FeQUAD4<typename T_MATRIX_TYPES::NodalVectorType, typename T_MATRIX_TYPES::DimNodalMatrixType, typename T_MATRIX_TYPES::DimMatrixType>;
    typedef MeshLib::Quad MeshElementType;
    static const unsigned dim = 2; //MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;

    MeshLib::Quad* createMeshElement()
    {
        // square
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0);
        return new MeshLib::Quad(nodes);
    }

    // set an expected mass matrix for 1m x 1m
    template <class T_MATRIX, typename ID_TYPE=signed>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.0; m(0,1) = 1./2; m(0,2) = 1./4; m(0,3) = 1./2;
        m(1,1) = 1.0; m(1,2) = 1./2; m(1,3) = 1./4;
        m(2,2) = 1.0; m(2,3) = 1./2;
        m(3,3) = 1.0;
        // make symmetric
        copyUpperToLower(4, m);
        m *= 1./9.;
    }

    // set an expected laplace matrix for 1m x 1m
    template <class T_MATRIX, typename ID_TYPE=signed>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 4.0; m(0,1) = -1.0; m(0,2) = -2.0; m(0,3) = -1.0;
        m(1,1) = 4.0; m(1,2) = -1.0; m(1,3) = -2.0;
        m(2,2) = 4.0; m(2,3) = -1.0;
        m(3,3) = 4.0;
        // make symmetric
        copyUpperToLower(4, m);
        m *= k/6.;
    }
};

// Test case for HEX8
class TestFeHEX8
{
public:
    template <class T_MATRIX_TYPES>
    using FeType = NumLib::FeHEX8<typename T_MATRIX_TYPES::NodalVectorType, typename T_MATRIX_TYPES::DimNodalMatrixType, typename T_MATRIX_TYPES::DimMatrixType>;
    typedef MeshLib::Hex MeshElementType;
    static const unsigned dim = 3; //MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;

    MeshLib::Hex* createMeshElement()
    {
        // cubic
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0);
        nodes[4] = new MeshLib::Node(0.0, 0.0, 1.0);
        nodes[5] = new MeshLib::Node(1.0, 0.0, 1.0);
        nodes[6] = new MeshLib::Node(1.0, 1.0, 1.0);
        nodes[7] = new MeshLib::Node(0.0, 1.0, 1.0);
        return new MeshLib::Hex(nodes);
    }

    // set an expected mass matrix
    template <class T_MATRIX, typename ID_TYPE=signed>
    void setExpectedMassMatrix(T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 1.0; m(0,1) = 1./2; m(0,2) = 1./4; m(0,3) = 1./2; m(0,4) = 1./2; m(0,5) = 1./4; m(0,6) = 1./8; m(0,7) = 1./4;
        m(1,1) = 1.0; m(1,2) = 1./2; m(1,3) = 1./4; m(1,4) = 1./4; m(1,5) = 1./2; m(1,6) = 1./4; m(1,7) = 1./8;
        m(2,2) = 1.0; m(2,3) = 1./2; m(2,4) = 1./8; m(2,5) = 1./4; m(2,6) = 1./2; m(2,7) = 1./4;
        m(3,3) = 1.0; m(3,4) = 1./4; m(3,5) = 1./8; m(3,6) = 1./4; m(3,7) = 1./2;
        m(4,4) = 1.0; m(4,5) = 1./2; m(4,6) = 1./4; m(4,7) = 1./2;
        m(5,5) = 1.0; m(5,6) = 1./2; m(5,7) = 1./4;
        m(6,6) = 1.0; m(6,7) = 1./2;
        m(7,7) = 1.0;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= 1./27.;
    }

    // set an expected laplace matrix
    template <class T_MATRIX, typename ID_TYPE=signed>
    void setExpectedLaplaceMatrix(double k, T_MATRIX &m)
    {
        // set upper triangle entries
        m(0,0) = 2./3; m(0,1) = 0; m(0,2) = -1./6; m(0,3) = 0; m(0,4) = 0; m(0,5) = -1./6; m(0,6) = -1./6; m(0,7) = -1./6;
        m(1,1) = 2./3; m(1,2) = 0; m(1,3) = -1./6; m(1,4) = -1./6; m(1,5) = 0; m(1,6) = -1./6; m(1,7) = -1./6;
        m(2,2) = 2./3; m(2,3) = 0; m(2,4) = -1./6; m(2,5) = -1./6; m(2,6) = 0; m(2,7) = -1./6;
        m(3,3) = 2./3; m(3,4) = -1./6; m(3,5) = -1./6; m(3,6) = -1./6; m(3,7) = 0;
        m(4,4) = 2./3; m(4,5) = 0; m(4,6) = -1./6; m(4,7) = 0;
        m(5,5) = 2./3; m(5,6) = 0; m(5,7) = -1./6;
        m(6,6) = 2./3; m(6,7) = 0;
        m(7,7) = 2./3;
        // make symmetric
        copyUpperToLower(e_nnodes, m);
        m *= k/2.;
    }
};

template <class T>
class NumLibFemIsoTest : public ::testing::Test, public T::T_FE
{
 public:
    typedef typename T::T_MATRIX_TYPES T_MATRIX_TYPES;
    typedef typename T::T_FE T_FE;
    // Matrix types
    typedef T_MATRIX_TYPES MatrixType;
    typedef typename T_MATRIX_TYPES::NodalMatrixType NodalMatrix;
    typedef typename T_MATRIX_TYPES::NodalVectorType NodalVector;
    typedef typename T_MATRIX_TYPES::DimNodalMatrixType DimNodalMatrix;
    typedef typename T_MATRIX_TYPES::DimMatrixType DimMatrix;
    // Finite element type
    typedef typename T_FE::template FeType<T_MATRIX_TYPES>::type FeType;
    // Shape matrix data type
    typedef typename FeType::ShapeMatricesType ShapeMatricesType;
    typedef typename T_FE::MeshElementType MeshElementType;

    static const unsigned dim = T_FE::dim;
    static const unsigned e_nnodes = T_FE::e_nnodes;

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
    typename FeType::IntegrationMethod integration_method;

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

} // namespace

#ifdef OGS_USE_EIGEN
template <typename T_FE>
struct EigenFixedMatrixTypes
{
    typedef Eigen::Matrix<double, T_FE::e_nnodes, T_FE::e_nnodes, Eigen::RowMajor> NodalMatrixType;
    typedef Eigen::Matrix<double, T_FE::e_nnodes, 1> NodalVectorType;
    typedef Eigen::Matrix<double, T_FE::dim, T_FE::e_nnodes, Eigen::RowMajor> DimNodalMatrixType;
    typedef Eigen::Matrix<double, T_FE::dim, T_FE::dim, Eigen::RowMajor> DimMatrixType;
};

struct EigenDynamicMatrixTypes
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> NodalMatrixType;
    typedef NodalMatrixType DimNodalMatrixType;
    typedef NodalMatrixType DimMatrixType;
    typedef Eigen::VectorXd NodalVectorType;
};

typedef EigenDynamicMatrixTypes DefaultMatrixType;
#endif // OGS_USE_EIGEN

// test cases
template <class T_FE_TYPE, class T_MAT_TYPES=DefaultMatrixType>
struct TestCase
{
    typedef T_FE_TYPE T_FE;
    typedef T_MAT_TYPES T_MATRIX_TYPES;
};

typedef ::testing::Types<
        TestCase<TestFeLINE2>,
        TestCase<TestFeQUAD4>,
        TestCase<TestFeHEX8>
#ifdef OGS_USE_EIGEN
        ,TestCase<TestFeLINE2, EigenFixedMatrixTypes<TestFeLINE2> >
        ,TestCase<TestFeQUAD4, EigenFixedMatrixTypes<TestFeQUAD4> >
        ,TestCase<TestFeHEX8, EigenFixedMatrixTypes<TestFeHEX8> >
#endif
        > TestTypes;


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
    NodalMatrix K(this->e_nnodes, this->e_nnodes);
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
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    ASSERT_EQ(std::pow(2u, this->dim), this->integration_method.getNPoints());
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
    ASSERT_EQ(std::pow(3u, this->dim), this->integration_method.getNPoints());
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.computeShapeFunctions(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}


