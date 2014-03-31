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

static const unsigned dim = 1;
static const unsigned e_nnodes = 2;

template <class T_MATRIX_TYPES>
class NumLibFemIsoLine2Test : public ::testing::Test
{
 public:
    // Matrix types
    typedef typename T_MATRIX_TYPES::NodalMatrixType NodalMatrix;
    typedef typename T_MATRIX_TYPES::NodalVectorType NodalVector;
    typedef typename T_MATRIX_TYPES::DimNodalMatrixType DimNodalMatrix;
    typedef typename T_MATRIX_TYPES::DimMatrixType DimMatrix;
    // Finite element type
    typedef typename NumLib::FeLINE2<NodalVector, DimNodalMatrix, DimMatrix>::type FeLINE2Type;
    // Shape matrix data type
    typedef typename FeLINE2Type::ShapeMatricesType ShapeMatricesType;

 public:
    NumLibFemIsoLine2Test() :
        D(dim, dim),
        expectedM(e_nnodes,e_nnodes),
        expectedK(e_nnodes,e_nnodes),
        integration_method(2)
    {
        // create a quad element used for testing
        unitLine = createLine(1.0);

        // set a conductivity tensor
        setIdentityMatrix(dim, D);
        D *= conductivity;

        // set expected matrices
        setExpectedMassMatrix(expectedM);
        setExpectedLaplaceMatrix(conductivity, expectedK);

        // for destructor
        vec_eles.push_back(unitLine);
        for (auto e : vec_eles)
            for (unsigned i=0; i<e->getNNodes(true); i++)
                vec_nodes.push_back(e->getNode(i));
    }

    virtual ~NumLibFemIsoLine2Test()
    {
        for (auto itr = vec_nodes.begin(); itr!=vec_nodes.end(); ++itr )
            delete *itr;
        for (auto itr = vec_eles.begin(); itr!=vec_eles.end(); ++itr )
            delete *itr;
    }

    // Line: a length of 1m
    MeshLib::Line* createLine(double h)
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(  h, 0.0, 0.0);
        return new MeshLib::Line(nodes);
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
        m(0,0) = 1./3.; m(0,1) = 1./6.;
        m(1,1) = 1./3.;
        // make symmetric
        copyUpperToLower(2, m);
    }

    // set an expected laplace matrix for 1m x 1m
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

    static const double conductivity;
    static const double eps;
    DimMatrix D;
    NodalMatrix expectedM;
    NodalMatrix expectedK;
    typename FeLINE2Type::IntegrationMethod integration_method;

    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const MeshLib::Line*> vec_eles;
    MeshLib::Line* unitLine;

}; // NumLibFemIsoLineTest

template <class T_MATRIX_TYPES>
const double NumLibFemIsoLine2Test<T_MATRIX_TYPES>::conductivity = 1e-11;

template <class T_MATRIX_TYPES>
const double NumLibFemIsoLine2Test<T_MATRIX_TYPES>::eps = std::numeric_limits<double>::epsilon();

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

TYPED_TEST_CASE(NumLibFemIsoLine2Test, MatrixTypes);

TYPED_TEST(NumLibFemIsoLine2Test, CheckMassMatrix)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeLINE2Type FeLINE2Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeLINE2Type fe(*this->unitLine);

    // evaluate a mass matrix M = int{ N^T D N }dA_e
    NodalMatrix M(e_nnodes, e_nnodes);
    ShapeMatricesType shape(dim, e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.template computeShapeFunctions<ShapeMatrixType::N_J>(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }

    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoLine2Test, CheckLaplaceMatrix)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeLINE2Type FeLINE2Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeLINE2Type fe(*this->unitLine);

    // evaluate a Laplace matrix K = int{ dNdx^T D dNdx }dA_e
    NodalMatrix K(e_nnodes, e_nnodes);
    ShapeMatricesType shape(dim, e_nnodes);
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.template computeShapeFunctions<ShapeMatrixType::DNDX>(wp.getCoords(), shape);
        K.noalias() += shape.dNdx.transpose() * this->D * shape.dNdx * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedK.data(), K.data(), K.size(), this->eps);
}

TYPED_TEST(NumLibFemIsoLine2Test, CheckMassLaplaceMatrices)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeLINE2Type FeLINE2Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object
    FeLINE2Type fe(*this->unitLine);

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

TYPED_TEST(NumLibFemIsoLine2Test, CheckGaussIntegrationLevel)
{
    // Refer to typedefs in the fixture
    typedef typename TestFixture::FeLINE2Type FeLINE2Type;
    typedef typename TestFixture::NodalMatrix NodalMatrix;
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;

    // create a finite element object with gauss quadrature level 2
    FeLINE2Type fe(*this->unitLine);

    // evaluate a mass matrix
    NodalMatrix M(e_nnodes, e_nnodes);
    ShapeMatricesType shape(dim, e_nnodes);
    ASSERT_EQ(2u, this->integration_method.getNPoints());
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
    ASSERT_EQ(3u, this->integration_method.getNPoints());
    for (std::size_t i=0; i < this->integration_method.getNPoints(); i++) {
        shape.setZero();
        auto wp = this->integration_method.getWeightedPoint(i);
        fe.computeShapeFunctions(wp.getCoords(), shape);
        M.noalias() += shape.N * shape.N.transpose() * shape.detJ * wp.getWeight();
    }
    ASSERT_ARRAY_NEAR(this->expectedM.data(), M.data(), M.size(), this->eps);
}


