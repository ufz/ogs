/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <limits>
#include <algorithm>
#include <vector>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Hex.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalCoordinatesMapping.h"

#include "../TestTools.h"

using namespace NumLib;

#ifdef OGS_USE_EIGEN

namespace
{

class Line2Test
{
public:
    typedef MeshLib::Line ElementType;
    typedef ShapeLine2 ShapeFunctionType;
    static const unsigned dim = ElementType::dimension;
    static const unsigned e_nnodes = ElementType::n_all_nodes;
    static const double r[dim];
    static const double nat_exp_N[e_nnodes];
    static const double nat_exp_dNdr[e_nnodes*dim];
    static const double ir_exp_J[dim*dim];
    static const double ir_exp_invJ[dim*dim];
    static const double ir_exp_detJ;
    static const double ir_exp_dNdx[e_nnodes*dim];
    static const double cl_exp_J[dim*dim];
    static const double cl_exp_detJ;
    static const double ze_exp_J[dim*dim];

    // element shape identical to that in natural coordinates (see ShapeLine2.h)
    static MeshLib::Line* createNaturalShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node( 1.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // element having irregular or skew shape
    MeshLib::Line* createIrregularShape()
    {
        // two times longer than the natural
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-2.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node( 2.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // invalid case: clock wise node ordering
    MeshLib::Line* createClockWise()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[1] = new MeshLib::Node(-1.0, 0.0, 0.0);
        nodes[0] = new MeshLib::Node( 1.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }

    // invalid case: zero volume
    MeshLib::Line* createZeroVolume()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(0.0, 0.0, 0.0);
        return new MeshLib::Line(nodes);
    }
};

const double Line2Test::r[dim] = {0.5};
const double Line2Test::nat_exp_N[e_nnodes] = {0.25, 0.75};
const double Line2Test::nat_exp_dNdr[e_nnodes*dim] = {-0.5, 0.5};
const double Line2Test::ir_exp_J[dim*dim]= {2.0};
const double Line2Test::ir_exp_invJ[dim*dim]= {0.5};
const double Line2Test::ir_exp_detJ = 2.;
const double Line2Test::ir_exp_dNdx[dim*e_nnodes]= {-0.25, 0.25};
const double Line2Test::cl_exp_J[dim*dim]= {-1.};
const double Line2Test::cl_exp_detJ = -1;
const double Line2Test::ze_exp_J[dim*dim]= {0.0};

class Quad4Test
{
 public:
    // Matrix types
    typedef MeshLib::Quad ElementType;
    typedef ShapeQuad4 ShapeFunctionType;
    static const unsigned dim = 2; //ElementType::dimension;
    static const unsigned e_nnodes = ElementType::n_all_nodes;

    // quad having shape identical to that in natural coordinates
    MeshLib::Quad* createNaturalShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 1.0,  1.0,  0.0);
        nodes[1] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[2] = new MeshLib::Node(-1.0, -1.0,  0.0);
        nodes[3] = new MeshLib::Node( 1.0, -1.0,  0.0);
        return new MeshLib::Quad(nodes);
    }

    // quad having irregular or skew shape
    MeshLib::Quad* createIrregularShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-0.5, -0.5,  0.0);
        nodes[1] = new MeshLib::Node( 0.6, -0.6,  0.0);
        nodes[2] = new MeshLib::Node( 0.5,  0.4,  0.0);
        nodes[3] = new MeshLib::Node(-0.3,  0.1,  0.0);
        return new MeshLib::Quad(nodes);
    }

    // invalid case: clock wise node ordering
    MeshLib::Quad* createClockWise()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 1.0,  1.0,  0.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[2] = new MeshLib::Node(-1.0, -1.0,  0.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0,  0.0);
        return new MeshLib::Quad(nodes);
    }

    // invalid case: zero area
    MeshLib::Quad* createZeroVolume()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node( 1.0,  1.0,  0.0);
        nodes[1] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[2] = new MeshLib::Node(-1.0,  1.0,  0.0);
        nodes[3] = new MeshLib::Node( 1.0,  1.0,  0.0);
        return new MeshLib::Quad(nodes);
    }

    static const double r[dim];
    static const double nat_exp_N[e_nnodes];
    static const double nat_exp_dNdr[e_nnodes*dim];
    static const double ir_exp_J[dim*dim];
    static const double ir_exp_invJ[dim*dim];
    static const double ir_exp_detJ;
    static const double ir_exp_dNdx[e_nnodes*dim];
    static const double cl_exp_J[dim*dim];
    static const double cl_exp_detJ;
    static const double ze_exp_J[dim*dim];
};

const double Quad4Test::r[dim] = {0.5, 0.5};
const double Quad4Test::nat_exp_N[e_nnodes] = {0.5625, 0.1875, 0.0625, 0.1875};
const double Quad4Test::nat_exp_dNdr[e_nnodes*dim] = {0.375, -0.375, -0.125, 0.125, 0.375, 0.125, -0.125, -0.375};
const double Quad4Test::ir_exp_J[dim*dim]= {-0.5125, 0.0, -0.0625, -0.35};
const double Quad4Test::ir_exp_detJ= 0.179375;
const double Quad4Test::ir_exp_invJ[dim*dim]= {-1.9512195121951219, 0.0, 0.3484320557491290, -2.8571428571428572};
const double Quad4Test::ir_exp_dNdx[dim*e_nnodes]= {-0.73170731707317072, 0.73170731707317072, 0.243902439024390, -0.24390243902439029, -0.940766550522648, -0.48780487804878048, 0.313588850174216, 1.1149825783972125};
const double Quad4Test::cl_exp_J[dim*dim]= {0.0, 1.0, 1.0, 0.0};
const double Quad4Test::cl_exp_detJ= -1.;
const double Quad4Test::ze_exp_J[dim*dim]= {1.0, 0.0, 0.0, 0.0};

class Hex8Test
{
 public:
    // Matrix types
    typedef MeshLib::Hex ElementType;
    typedef ShapeHex8 ShapeFunctionType;
    static const unsigned dim = 3; //ElementType::dimension;
    static const unsigned e_nnodes = ElementType::n_all_nodes;

    // element shape identical to that in natural coordinates (see ShapeHex8.h)
    MeshLib::Hex* createNaturalShape()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0, -1.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0, -1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0, -1.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0, -1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[5] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[7] = new MeshLib::Node(-1.0,  1.0,  1.0);
        return new MeshLib::Hex(nodes);
    }

    // element having irregular or skew shape
    MeshLib::Hex* createIrregularShape()
    {
        // two times longer in z direction than the natural
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0, -1.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0, -1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0, -1.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0, -1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  3.0);
        nodes[5] = new MeshLib::Node( 1.0, -1.0,  3.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  3.0);
        nodes[7] = new MeshLib::Node(-1.0,  1.0,  3.0);
        return new MeshLib::Hex(nodes);
    }

    // invalid case: clock wise node ordering
    MeshLib::Hex* createClockWise()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0, -1.0);
        nodes[3] = new MeshLib::Node( 1.0, -1.0, -1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0, -1.0);
        nodes[1] = new MeshLib::Node(-1.0,  1.0, -1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[7] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[5] = new MeshLib::Node(-1.0,  1.0,  1.0);
        return new MeshLib::Hex(nodes);
    }

    // invalid case: zero volume
    MeshLib::Hex* createZeroVolume()
    {
        MeshLib::Node** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[1] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[2] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[3] = new MeshLib::Node(-1.0,  1.0,  1.0);
        nodes[4] = new MeshLib::Node(-1.0, -1.0,  1.0);
        nodes[5] = new MeshLib::Node( 1.0, -1.0,  1.0);
        nodes[6] = new MeshLib::Node( 1.0,  1.0,  1.0);
        nodes[7] = new MeshLib::Node(-1.0,  1.0,  1.0);
        return new MeshLib::Hex(nodes);
    }

    static const double r[dim];
    static const double nat_exp_N[e_nnodes];
    static const double nat_exp_dNdr[e_nnodes*dim];
    static const double ir_exp_J[dim*dim];
    static const double ir_exp_invJ[dim*dim];
    static const double ir_exp_detJ;
    static const double ir_exp_dNdx[e_nnodes*dim];
    static const double cl_exp_J[dim*dim];
    static const double cl_exp_detJ;
    static const double ze_exp_J[dim*dim];
};

const double Hex8Test::r[dim] = {0.5, 0.5, 0.5};
const double Hex8Test::nat_exp_N[e_nnodes]
    = {0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625};
const double Hex8Test::nat_exp_dNdr[e_nnodes*dim]
    = {-0.03125, 0.03125, 0.09375, -0.09375, -0.09375, 0.09375, 0.28125, -0.28125,
       -0.03125, -0.09375, 0.09375, 0.03125, -0.09375, -0.28125,  0.28125,  0.09375,
       -0.03125, -0.09375, -0.28125, -0.09375,  0.03125,  0.09375,  0.28125,  0.09375};
const double Hex8Test::ir_exp_J[dim*dim]= {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0};
const double Hex8Test::ir_exp_detJ= 2.;
const double Hex8Test::ir_exp_invJ[dim*dim]= {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1./2.0};
const double Hex8Test::ir_exp_dNdx[dim*e_nnodes]
    = {-0.03125, 0.03125, 0.09375, -0.09375, -0.09375, 0.09375, 0.28125, -0.28125,
       -0.03125, -0.09375, 0.09375, 0.03125, -0.09375, -0.28125,  0.28125,  0.09375,
       -0.015625, -0.046875, -0.140625, -0.046875,  0.015625,  0.046875,  0.140625,  0.046875};
const double Hex8Test::cl_exp_J[dim*dim]= {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
const double Hex8Test::cl_exp_detJ= -1.;
const double Hex8Test::ze_exp_J[dim*dim]= {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};


template <class T_TEST>
class NumLibFemNaturalCoordinatesMappingTest : public ::testing::Test, public T_TEST
{
public:
    typedef typename T_TEST::ElementType ElementType;
    typedef typename T_TEST::ShapeFunctionType ShapeFunctionType;
    static const unsigned dim = T_TEST::dim;
    static const unsigned e_nnodes = T_TEST::e_nnodes;
    // Matrix types
    typedef Eigen::Matrix<double, e_nnodes, 1> NodalVector;
    typedef Eigen::Matrix<double, dim, e_nnodes, Eigen::RowMajor> DimNodalMatrix;
    typedef Eigen::Matrix<double, dim, dim, Eigen::RowMajor> DimMatrix;
    // Shape data type
    typedef ShapeMatrices<NodalVector,DimNodalMatrix,DimMatrix> ShapeMatricesType;
    // Natural coordinates mapping type
    typedef NaturalCoordinatesMapping<ElementType, ShapeFunctionType, ShapeMatricesType> NaturalCoordsMappingType;

public:
    NumLibFemNaturalCoordinatesMappingTest() : eps(std::numeric_limits<double>::epsilon())
    {
        // create four elements used for testing
        naturalEle   = this->createNaturalShape();
        irregularEle = this->createIrregularShape();
        clockwiseEle  = this->createClockWise();
        zeroVolumeEle  = this->createZeroVolume();

        // for destructor
        vec_eles.push_back(naturalEle);
        vec_eles.push_back(irregularEle);
        vec_eles.push_back(clockwiseEle);
        vec_eles.push_back(zeroVolumeEle);
        for (auto e : vec_eles)
            for (unsigned i=0; i<e->getNNodes(true); i++)
                vec_nodes.push_back(e->getNode(i));
    }

    ~NumLibFemNaturalCoordinatesMappingTest()
    {
        for (auto itr = vec_nodes.begin(); itr!=vec_nodes.end(); ++itr )
            delete *itr;
        for (auto itr = vec_eles.begin(); itr!=vec_eles.end(); ++itr )
            delete *itr;
    }

    const double eps;
    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const ElementType*> vec_eles;
    ElementType* naturalEle;
    ElementType* irregularEle;
    ElementType* clockwiseEle;
    ElementType* zeroVolumeEle;
};

} // namespace

typedef ::testing::Types<
        Line2Test,
        Quad4Test,
        Hex8Test
        > ElementTypes;

TYPED_TEST_CASE(NumLibFemNaturalCoordinatesMappingTest, ElementTypes);

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_N)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->e_nnodes);

    //only N
    NaturalCoordsMappingType::template computeShapeMatrices<ShapeMatrixType::N>(*this->naturalEle, this->r, shape);
    ASSERT_FALSE(shape.N.isZero());
    ASSERT_TRUE(shape.dNdr.isZero());
    ASSERT_TRUE(shape.J.isZero());
    ASSERT_TRUE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_DNDR)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->e_nnodes);

    // dNdr
    NaturalCoordsMappingType::template computeShapeMatrices<ShapeMatrixType::DNDR>(*this->naturalEle, this->r, shape);
    ASSERT_TRUE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_TRUE(shape.J.isZero());
    ASSERT_TRUE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_N_J)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->e_nnodes);

    // N_J
    shape.setZero();
    NaturalCoordsMappingType::template computeShapeMatrices<ShapeMatrixType::N_J>(*this->naturalEle, this->r, shape);
    ASSERT_FALSE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_DNDR_J)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->e_nnodes);

    // dNdr, J
    NaturalCoordsMappingType::template computeShapeMatrices<ShapeMatrixType::DNDR_J>(*this->naturalEle, this->r, shape);
    ASSERT_TRUE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.isZero());
    ASSERT_TRUE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_DNDX)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->e_nnodes);

    // DNDX
    shape.setZero();
    NaturalCoordsMappingType::template computeShapeMatrices<ShapeMatrixType::DNDX>(*this->naturalEle, this->r, shape);
    ASSERT_TRUE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_FALSE(shape.invJ.isZero());
    ASSERT_FALSE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckFieldSpecification_ALL)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;
    ShapeMatricesType shape(this->dim, this->e_nnodes);

    // ALL
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices(*this->naturalEle, this->r, shape);
    ASSERT_FALSE(shape.N.isZero());
    ASSERT_FALSE(shape.dNdr.isZero());
    ASSERT_FALSE(shape.J.isZero());
    ASSERT_FALSE(shape.detJ == .0);
    ASSERT_FALSE(shape.invJ.isZero());
    ASSERT_FALSE(shape.dNdx.isZero());
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckNaturalShape)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;

    // identical to natural coordinates
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    NaturalCoordsMappingType::computeShapeMatrices(*this->naturalEle, this->r, shape);
    double exp_J[this->dim*this->dim]= {0.0};
    for (unsigned i=0; i<this->dim; i++)
        exp_J[i+this->dim*i] = 1.0;

    ASSERT_ARRAY_NEAR(this->nat_exp_N, shape.N.data(), shape.N.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), this->eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), this->eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.invJ.data(), shape.invJ.size(), this->eps);
    ASSERT_NEAR(1.0, shape.detJ, this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdx.data(), shape.dNdx.size(), this->eps);
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckIrregularShape)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;

    // irregular shape
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    NaturalCoordsMappingType::computeShapeMatrices(*this->irregularEle, this->r, shape);
    //std::cout << shape;

    ASSERT_ARRAY_NEAR(this->nat_exp_N, shape.N.data(), shape.N.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->ir_exp_J, shape.J.data(), shape.J.size(), this->eps);
    ASSERT_NEAR(this->ir_exp_detJ, shape.detJ, this->eps);
    ASSERT_ARRAY_NEAR(this->ir_exp_invJ, shape.invJ.data(), shape.invJ.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->ir_exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), this->eps);
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckClockwise)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;

    // clockwise node ordering, which is invalid)
    ShapeMatricesType shape(this->dim, this->e_nnodes);
    NaturalCoordsMappingType::computeShapeMatrices(*this->clockwiseEle, this->r, shape);
    //std::cout << shape;
    // Inverse of the Jacobian matrix doesn't exist
    double exp_invJ[this->dim*this->dim]= {0.0};
    double exp_dNdx[this->dim*this->e_nnodes]= {0.0};

    ASSERT_ARRAY_NEAR(this->nat_exp_N, shape.N.data(), shape.N.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->cl_exp_J, shape.J.data(), shape.J.size(), this->eps);
    ASSERT_NEAR(this->cl_exp_detJ, shape.detJ, this->eps);
    ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), this->eps);
    ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), this->eps);
}

TYPED_TEST(NumLibFemNaturalCoordinatesMappingTest, CheckZeroVolume)
{
    typedef typename TestFixture::ShapeMatricesType ShapeMatricesType;
    typedef typename TestFixture::NaturalCoordsMappingType NaturalCoordsMappingType;

    ShapeMatricesType shape(this->dim, this->e_nnodes);
    NaturalCoordsMappingType::computeShapeMatrices(*this->zeroVolumeEle, this->r, shape);
    //std::cout << shape;
    // Inverse of the Jacobian matrix doesn't exist
    double exp_invJ[this->dim*this->dim]= {0.0};
    double exp_dNdx[this->dim*this->e_nnodes]= {0.0};

    ASSERT_ARRAY_NEAR(this->nat_exp_N, shape.N.data(), shape.N.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->nat_exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), this->eps);
    ASSERT_ARRAY_NEAR(this->ze_exp_J, shape.J.data(), shape.J.size(), this->eps);
    ASSERT_NEAR(0.0, shape.detJ, this->eps);
    ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), this->eps);
    ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), this->eps);
}

#endif //OGS_USE_EIGEN


