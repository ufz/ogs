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
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalCoordinatesMapping.h"

#include "../TestTools.h"

using namespace NumLib;

#ifdef OGS_USE_EIGEN

namespace
{
class NumLibFemNaturalCoordinatesMappingLine2Test : public ::testing::Test
{
 public:
    // Matrix types
    static const unsigned dim = 1;
    static const unsigned e_nnodes = 2;
    typedef Eigen::Matrix<double, e_nnodes, 1> NodalVector;
    typedef Eigen::Matrix<double, dim, e_nnodes, Eigen::RowMajor> DimNodalMatrix;
    typedef Eigen::Matrix<double, dim, dim, Eigen::RowMajor> DimMatrix;
    // Shape data type
    typedef ShapeMatrices<NodalVector,DimNodalMatrix,DimMatrix> ShapeMatricesType;
    // Natural coordinates mapping type
    typedef NaturalCoordinatesMapping<MeshLib::Line, ShapeLine2, ShapeMatricesType> NaturalCoordsMappingType;

 public:
    NumLibFemNaturalCoordinatesMappingLine2Test()
    {
        // create four elements used for testing
        naturalEle   = createNaturalShape();
        irregularEle = createIrregularShape();
        clockwiseEle  = createClockWise();
        zeroVolumeEle  = createZeroVolume();

        // for destructor
        vec_eles.push_back(naturalEle);
        vec_eles.push_back(irregularEle);
        vec_eles.push_back(clockwiseEle);
        vec_eles.push_back(zeroVolumeEle);
        for (auto e : vec_eles)
            for (unsigned i=0; i<e->getNNodes(true); i++)
                vec_nodes.push_back(e->getNode(i));
    }

    ~NumLibFemNaturalCoordinatesMappingLine2Test()
    {
        for (auto itr = vec_nodes.begin(); itr!=vec_nodes.end(); ++itr )
            delete *itr;
        for (auto itr = vec_eles.begin(); itr!=vec_eles.end(); ++itr )
            delete *itr;
    }

    // element shape identical to that in natural coordinates (see ShapeLine2.h)
    MeshLib::Line* createNaturalShape()
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

    static const double r[dim];
    static const double exp_N[e_nnodes];
    static const double exp_dNdr[e_nnodes*dim];
    static const double eps;

    std::vector<const MeshLib::Node*> vec_nodes;
    std::vector<const MeshLib::Line*> vec_eles;
    MeshLib::Line* naturalEle;
    MeshLib::Line* irregularEle;
    MeshLib::Line* clockwiseEle;
    MeshLib::Line* zeroVolumeEle;

}; // NumLibFemNaturalCoordinatesMappingHex8Test

const double NumLibFemNaturalCoordinatesMappingLine2Test::r[dim] = {0.5};
const double NumLibFemNaturalCoordinatesMappingLine2Test::exp_N[e_nnodes] = {0.25, 0.75};
const double NumLibFemNaturalCoordinatesMappingLine2Test::exp_dNdr[e_nnodes*dim] = {-0.5, 0.5};
const double NumLibFemNaturalCoordinatesMappingLine2Test::eps = std::numeric_limits<double>::epsilon();

} // namespace

TEST_F(NumLibFemNaturalCoordinatesMappingLine2Test, CheckNaturalShape)
{
    // identical to natural coordinates
    ShapeMatricesType shape(dim, e_nnodes);
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices(*naturalEle, r, shape);
    double exp_J[]= {1.0};

    ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.invJ.data(), shape.invJ.size(), eps);
    ASSERT_NEAR(1.0, shape.detJ, eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdx.data(), shape.dNdx.size(), eps);
}

TEST_F(NumLibFemNaturalCoordinatesMappingLine2Test, CheckIrregularShape)
{
    // irregular shape
    ShapeMatricesType shape(dim, e_nnodes);
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices(*irregularEle, r, shape);
    //std::cout << shape;
    double exp_J[]= {2.0};
    double exp_invJ[]= {1./2.0};
    double exp_dNdx[] = {-0.25, 0.25};

    ASSERT_ARRAY_NEAR(exp_N, shape.N.data(), shape.N.size(), eps);
    ASSERT_ARRAY_NEAR(exp_dNdr, shape.dNdr.data(), shape.dNdr.size(), eps);
    ASSERT_ARRAY_NEAR(exp_J, shape.J.data(), shape.J.size(), eps);
    ASSERT_ARRAY_NEAR(exp_invJ, shape.invJ.data(), shape.invJ.size(), eps);
    ASSERT_NEAR(2.0, shape.detJ, eps);
    ASSERT_ARRAY_NEAR(exp_dNdx, shape.dNdx.data(), shape.dNdx.size(), eps);
}

TEST_F(NumLibFemNaturalCoordinatesMappingLine2Test, CheckClockwise)
{
    // clockwise node ordering, which is invalid)
    ShapeMatricesType shape(dim, e_nnodes);
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices(*clockwiseEle, r, shape);
//    std::cout << shape;
    double exp_J[]= {-1.0};
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

TEST_F(NumLibFemNaturalCoordinatesMappingLine2Test, CheckZeroVolume)
{
    ShapeMatricesType shape(dim, e_nnodes);
    shape.setZero();
    NaturalCoordsMappingType::computeShapeMatrices(*zeroVolumeEle, r, shape);
//    std::cout << shape;
    double exp_J[dim*dim]= {0.0};
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


