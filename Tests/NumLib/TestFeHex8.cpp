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

#include "MeshLib/Elements/Hex.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

#include "Tests/TestTools.h"

using namespace NumLib;

namespace
{

static const unsigned dim = 3;
static const unsigned e_nnodes = 8;

template <class T_MATRIX_TYPES>
class NumLibFemIsoHex8Test : public ::testing::Test
{
 public:
    // Matrix types
    typedef typename T_MATRIX_TYPES::NodalMatrixType NodalMatrix;
    typedef typename T_MATRIX_TYPES::NodalVectorType NodalVector;
    typedef typename T_MATRIX_TYPES::DimNodalMatrixType DimNodalMatrix;
    typedef typename T_MATRIX_TYPES::DimMatrixType DimMatrix;
    // Finite element type
    typedef typename NumLib::FeHEX8<NodalVector, DimNodalMatrix, DimMatrix>::type FeHEX8Type;
    // Shape matrix data type
    typedef typename FeHEX8Type::ShapeMatricesType ShapeMatricesType;

 public:
    NumLibFemIsoHex8Test() :
        D(dim, dim),
        expectedM(e_nnodes,e_nnodes),
        expectedK(e_nnodes,e_nnodes),
        integration_method(2)
    {
        // create a quad element used for testing
        mesh_element   = createCubicHex(1.0);

        // set a conductivity tensor
        setIdentityMatrix(dim, D);
        D *= conductivity;

        // set expected matrices
        setExpectedMassMatrix(expectedM);
        setExpectedLaplaceMatrix(conductivity, expectedK);

        // for destructor
        vec_eles.push_back(mesh_element);
        for (auto e : vec_eles)
            for (unsigned i=0; i<e->getNNodes(true); i++)
                vec_nodes.push_back(e->getNode(i));
    }

    ~NumLibFemIsoHex8Test()
    {
        for (auto itr = vec_nodes.begin(); itr!=vec_nodes.end(); ++itr )
            delete *itr;
        for (auto itr = vec_eles.begin(); itr!=vec_eles.end(); ++itr )
            delete *itr;
    }

    // Quad: square shape with a length of 1m
    MeshLib::Hex* createCubicHex(double h)
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(  h, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(  h,   h, 0.0);
        nodes[3] = new MeshLib::Node(0.0,   h, 0.0);
        nodes[4] = new MeshLib::Node(0.0, 0.0, h);
        nodes[5] = new MeshLib::Node(  h, 0.0, h);
        nodes[6] = new MeshLib::Node(  h,   h, h);
        nodes[7] = new MeshLib::Node(0.0,   h, h);
        return new MeshLib::Hex(nodes);
    }

    // set an identity matrix
    template <class T_MATRIX, typename ID_TYPE=signed>
    void setIdentityMatrix(unsigned dim, T_MATRIX &m) const
    {
        for (unsigned i=0; i<dim; i++)
            for (unsigned j=0; j<dim; j++)
                m(i,j) = 0.0;
        for (unsigned i=0; i<dim; i++)
            m(i,i) = 1.0;
    }

    // copy upper triangles to lower triangles in a matrix
    template <class T_MATRIX, typename ID_TYPE=signed>
    void copyUpperToLower(const ID_TYPE dim, T_MATRIX &m) const
    {
        for (ID_TYPE i=0; i<dim; i++)
            for (ID_TYPE j=0; j<i; j++)
                m(i,j) = m(j,i);
    }

    // set an expected mass matrix for 1m x 1m
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

    // set an expected laplace matrix for 1m x 1m
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

    static const double conductivity;
    static const double eps;
    DimMatrix D;
    NodalMatrix expectedM;
    NodalMatrix expectedK;
    typename FeHEX8Type::IntegrationMethod integration_method;

    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const MeshLib::Hex*> vec_eles;
    MeshLib::Hex* mesh_element;

}; // NumLibFemIsoQuad4Test

template <class T_MATRIX_TYPES>
const double NumLibFemIsoHex8Test<T_MATRIX_TYPES>::conductivity = 1e-11;

template <class T_MATRIX_TYPES>
const double NumLibFemIsoHex8Test<T_MATRIX_TYPES>::eps = std::numeric_limits<double>::epsilon();

} // namespace

#ifdef OGS_USE_EIGEN

struct EigenFixedMatrixTypes
{
    typedef Eigen::Matrix<double, e_nnodes, e_nnodes, Eigen::RowMajor> NodalMatrixType;
    typedef Eigen::Matrix<double, e_nnodes, 1> NodalVectorType;
    typedef Eigen::Matrix<double, dim, e_nnodes, Eigen::RowMajor> DimNodalMatrixType;
    typedef Eigen::Matrix<double, dim, dim, Eigen::RowMajor> DimMatrixType;
};

struct EigenDynamicMatrixTypes
{
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> NodalMatrixType;
    typedef NodalMatrixType DimNodalMatrixType;
    typedef NodalMatrixType DimMatrixType;
    typedef Eigen::VectorXd NodalVectorType;
};

#endif // OGS_USE_EIGEN

typedef ::testing::Types<
#ifdef OGS_USE_EIGEN
        EigenFixedMatrixTypes
        , EigenDynamicMatrixTypes
#endif
        > MatrixTypes;

TYPED_TEST_CASE(NumLibFemIsoHex8Test, MatrixTypes);

TYPED_TEST(NumLibFemIsoHex8Test, CheckMassMatrix)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeHEX8Type FeHEX8Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeHEX8Type fe(*this->mesh_element);

    // evaluate a mass matrix M = int{ N^T D N }dA_e
    NodalMatrix M(e_nnodes, e_nnodes);
    ShapeMatricesType shape(dim, e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.template computeShapeFunctions<ShapeMatrixType::N_J>(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }

//    std::cout << this->expectedM << "\n M=";
//    std::cout << M;

    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoHex8Test, CheckLaplaceMatrix)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeHEX8Type FeHEX8Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeHEX8Type fe(*this->mesh_element);

    // evaluate a Laplace matrix K = int{ dNdx^T D dNdx }dA_e
    NodalMatrix K(e_nnodes, e_nnodes);
    ShapeMatricesType shape(dim, e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.template computeShapeFunctions<ShapeMatrixType::DNDX>(wp.getCoords(), shape);
        K.noalias() += shape.dNdx.transpose() * this->D * shape.dNdx * shape.detJ * wp.getWeight();
    }
//    std::cout << this->expectedK << "\n M=";
//    std::cout << K;
    ASSERT_ARRAY_NEAR(this->expectedK.data(), K.data(), K.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoHex8Test, CheckMassLaplaceMatrices)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeHEX8Type FeHEX8Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeHEX8Type fe(*this->mesh_element);

    // evaluate both mass and laplace matrices at once
    NodalMatrix M(e_nnodes, e_nnodes);
    NodalMatrix K(e_nnodes, e_nnodes);
    ShapeMatricesType shape(dim, e_nnodes);
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

TYPED_TEST(NumLibFemIsoHex8Test, CheckGaussIntegrationLevel)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeHEX8Type FeHEX8Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object with gauss quadrature level 2
    FeHEX8Type fe(*this->mesh_element);

    // evaluate a mass matrix
    NodalMatrix M(e_nnodes, e_nnodes);
    ShapeMatricesType shape(dim, e_nnodes);
    ASSERT_EQ(8u, this->integration_method.getNPoints());
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
    ASSERT_EQ(27u, this->integration_method.getNPoints());
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.computeShapeFunctions(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}


