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

#include <limits>
#include <algorithm>
#include <vector>

#ifdef OGS_USE_EIGEN
#include <Eigen>
#endif

#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalCoordinatesMapping.h"

#include "../TestTools.h"

using namespace NumLib;

#ifdef OGS_USE_EIGEN

namespace
{
class NumLibFemNaturalCoordinatesMappingQuad4Test : public ::testing::Test
{
 public:
    // Matrix types
    static const unsigned dim = 2;
    static const unsigned e_nnodes = 4;
    typedef Eigen::Matrix<double, e_nnodes, 1> NodalVector;
    typedef Eigen::Matrix<double, dim, e_nnodes, Eigen::RowMajor> DimNodalMatrix;
    typedef Eigen::Matrix<double, dim, dim, Eigen::RowMajor> DimMatrix;
    // Shape data type
    typedef ShapeMatrices<NodalVector,DimNodalMatrix,DimMatrix> ShapeMatricesType;
    // Natural coordinates mapping type
    typedef NaturalCoordinatesMapping<MeshLib::Quad, ShapeQuad4, ShapeMatricesType> NaturalCoordsMappingType;

 public:
    NumLibFemNaturalCoordinatesMappingQuad4Test()
    {
        // 0~3 natural
        vec_nodes.push_back(new MeshLib::Node( 1.0,  1.0,  0.0));
        vec_nodes.push_back(new MeshLib::Node(-1.0,  1.0,  0.0));
        vec_nodes.push_back(new MeshLib::Node(-1.0, -1.0,  0.0));
        vec_nodes.push_back(new MeshLib::Node( 1.0, -1.0,  0.0));
        // 4~7 irregular
        vec_nodes.push_back(new MeshLib::Node(-0.5, -0.5,  0.0));
        vec_nodes.push_back(new MeshLib::Node( 0.6, -0.6,  0.0));
        vec_nodes.push_back(new MeshLib::Node( 0.5,  0.4,  0.0));
        vec_nodes.push_back(new MeshLib::Node(-0.3,  0.1,  0.0));
        // 8~9 duplicated of node 0 and 1
        vec_nodes.push_back(new MeshLib::Node( 1.0,  1.0,  0.0));
        vec_nodes.push_back(new MeshLib::Node(-1.0,  1.0,  0.0));

        // quad having shape identical to that in natural coordiantes
        MeshLib::Node** naturalNodes = new MeshLib::Node*[e_nnodes];
        std::copy_n(vec_nodes.begin(), e_nnodes, naturalNodes);
        naturalQuad = new MeshLib::Quad(naturalNodes);
        vec_eles.push_back(naturalQuad);

        // irregular shape
        MeshLib::Node** irregularNodes = new MeshLib::Node*[e_nnodes];
        std::copy_n(vec_nodes.begin()+4, e_nnodes, irregularNodes);
        irregularQuad = new MeshLib::Quad(irregularNodes);
        vec_eles.push_back(irregularQuad);

        // invalid case: clock wise node ordering
        MeshLib::Node** clockwiseNodes = new MeshLib::Node*[e_nnodes];
        std::copy_n(vec_nodes.begin(), e_nnodes, naturalNodes);
        clockwiseNodes[0] = vec_nodes[0];
        clockwiseNodes[1] = vec_nodes[3];
        clockwiseNodes[2] = vec_nodes[2];
        clockwiseNodes[3] = vec_nodes[1];
        clockwiseQuad = new MeshLib::Quad(clockwiseNodes);
        vec_eles.push_back(clockwiseQuad);

        // invalid case: zero area
        MeshLib::Node** zeroAreaNodes = new MeshLib::Node*[e_nnodes];
        std::copy_n(vec_nodes.begin(), e_nnodes, naturalNodes);
        zeroAreaNodes[0] = vec_nodes[0];
        zeroAreaNodes[1] = vec_nodes[1];
        zeroAreaNodes[2] = vec_nodes[9];
        zeroAreaNodes[3] = vec_nodes[8];
        zeroAreaQuad = new MeshLib::Quad(zeroAreaNodes);
        vec_eles.push_back(zeroAreaQuad);
    }

    ~NumLibFemNaturalCoordinatesMappingQuad4Test()
    {
        for (auto itr = vec_nodes.begin(); itr!=vec_nodes.end(); ++itr )
            delete *itr;
        for (auto itr = vec_eles.begin(); itr!=vec_eles.end(); ++itr )
            delete *itr;
    }

    static const double r[dim];
    static const double exp_N[e_nnodes];
    static const double exp_dNdr[e_nnodes*dim];
    static const double eps;

    std::vector<MeshLib::Node*> vec_nodes;
    std::vector<MeshLib::Quad*> vec_eles;
    MeshLib::Quad* naturalQuad;
    MeshLib::Quad* irregularQuad;
    MeshLib::Quad* clockwiseQuad;
    MeshLib::Quad* zeroAreaQuad;

}; // NumLibFemNaturalCoordinatesMappingQuad4Test

const double NumLibFemNaturalCoordinatesMappingQuad4Test::r[dim] = {0.5, 0.5};
const double NumLibFemNaturalCoordinatesMappingQuad4Test::exp_N[e_nnodes] = {0.5625, 0.1875, 0.0625, 0.1875};
const double NumLibFemNaturalCoordinatesMappingQuad4Test::exp_dNdr[e_nnodes*dim] = {0.375, -0.375, -0.125, 0.125, 0.375, 0.125, -0.125, -0.375};
const double NumLibFemNaturalCoordinatesMappingQuad4Test::eps = std::numeric_limits<double>::epsilon();

} // namespace

TEST_F(NumLibFemNaturalCoordinatesMappingQuad4Test, CheckFieldSpecification)
{
    // identical to natural coordinates
    ShapeMatricesType shape(dim, e_nnodes);

    //only N
    NaturalCoordsMappingType::computeShapeMatrices<ShapeMatrixType::N>(*naturalQuad, r, shape);
    ASSERT_TRUE(shape.N.norm() > .0);
    ASSERT_TRUE(shape.dNdr.norm() == .0);
    ASSERT_TRUE(shape.J.norm() == .0);
    ASSERT_TRUE(shape.detJ == .0);
    ASSERT_TRUE(shape.invJ.norm() == .0);
    ASSERT_TRUE(shape.dNdx.norm() == .0);

    // N_J
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices<ShapeMatrixType::N_J>(*naturalQuad, r, shape);
    ASSERT_TRUE(shape.N.norm() > .0);
    ASSERT_TRUE(shape.dNdr.norm() > .0);
    ASSERT_TRUE(shape.J.norm() > .0);
    ASSERT_TRUE(shape.detJ > .0);
    ASSERT_TRUE(shape.invJ.norm() == .0);
    ASSERT_TRUE(shape.dNdx.norm() == .0);

    // DNDX
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices<ShapeMatrixType::DNDX>(*naturalQuad, r, shape);
    ASSERT_TRUE(shape.N.norm() == .0);
    ASSERT_TRUE(shape.dNdr.norm() > .0);
    ASSERT_TRUE(shape.J.norm() > .0);
    ASSERT_TRUE(shape.detJ > .0);
    ASSERT_TRUE(shape.invJ.norm() > .0);
    ASSERT_TRUE(shape.dNdx.norm() > .0);

    // ALL
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices(*naturalQuad, r, shape);
    ASSERT_TRUE(shape.N.norm() > .0);
    ASSERT_TRUE(shape.dNdr.norm() > .0);
    ASSERT_TRUE(shape.J.norm() > .0);
    ASSERT_TRUE(shape.detJ > .0);
    ASSERT_TRUE(shape.invJ.norm() > .0);
    ASSERT_TRUE(shape.dNdx.norm() > .0);
}


TEST_F(NumLibFemNaturalCoordinatesMappingQuad4Test, CheckNaturalShape)
{
    // identical to natural coordinates
    ShapeMatricesType shape(dim, e_nnodes);

    NaturalCoordsMappingType::computeShapeMatrices(*naturalQuad, r, shape);
    double exp_J[]= {1.0, 0.0, 0.0, 1.0};

    ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.invJ.data(), shape.invJ.size(), eps);
    ASSERT_NEAR(1.0, shape.detJ, eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdx.data(), shape.dNdx.size(), eps);
}

TEST_F(NumLibFemNaturalCoordinatesMappingQuad4Test, CheckIrregularShape)
{
    // irregular shape
    ShapeMatricesType shape(dim, e_nnodes);

    NaturalCoordsMappingType::computeShapeMatrices(*irregularQuad, r, shape);
//        std::cout << shape;
    double exp_J[]= {-0.5125, 0.0, -0.0625, -0.35};
    double exp_invJ[]= {-1.9512195121951219, 0.0, 0.3484320557491290, -2.8571428571428572};
    double exp_dNdx[]= {-0.73170731707317072, 0.73170731707317072, 0.243902439024390, -0.24390243902439029, -0.940766550522648, -0.48780487804878048, 0.313588850174216, 1.1149825783972125};

    ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
    ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), eps);
    ASSERT_NEAR(0.179375, shape.detJ, eps);
    ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), eps);
}

TEST_F(NumLibFemNaturalCoordinatesMappingQuad4Test, CheckClockwise)
{
    // clockwise node ordering, which is invalid)
    ShapeMatricesType shape(dim, e_nnodes);

    NaturalCoordsMappingType::computeShapeMatrices(*clockwiseQuad, r, shape);
    //std::cout << shape;
    double exp_J[]= {0.0, 1.0, 1.0, 0.0};
    // Inverse of the Jacobian matrix doesn't exist
    double exp_invJ[dim*dim]= {0.0};
    double exp_dNdx[dim*e_nnodes]= {0.0};

    ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
    ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), eps);
    ASSERT_NEAR(-1.0, shape.detJ, eps);
    ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), eps);
}

TEST_F(NumLibFemNaturalCoordinatesMappingQuad4Test, CheckZeroArea)
{
    // zero area
    ShapeMatricesType shape(dim, e_nnodes);

    NaturalCoordsMappingType::computeShapeMatrices(*zeroAreaQuad, r, shape);
    //std::cout << shape;
    double exp_J[]= {1.0, 0.0, 0.0, 0.0};
    // Inverse of the Jacobian matrix doesn't exist
    double exp_invJ[dim*dim]= {0.0};
    double exp_dNdx[dim*e_nnodes]= {0.0};

    ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
    ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), eps);
    ASSERT_NEAR(0.0, shape.detJ, eps);
    ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), eps);
}

#endif //OGS_USE_EIGEN


