/**
 * @brief Tests for GeoLib::Surface::isPntInSfc()
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <array>
#include <memory>
#include <random>
#include <vector>

#include "gtest/gtest.h"

#include <boost/math/constants/constants.hpp>

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"
#include "GeoLib/AnalyticalGeometry.h"

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/Point3d.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/convertMeshToGeo.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

inline double constant(double , double )
{
    return 0.0;
}

inline double coscos(double x, double y)
{
    return std::cos(x) * std::cos(y);
}

inline MathLib::Point3d
getEdgeMiddlePoint(MathLib::Point3d const&a, MathLib::Point3d const& b)
{
    return MathLib::Point3d(std::array<double,3>({{
        (a[0]+b[0])/2, (a[1]+b[1])/2, (a[2]+b[2])/2}}));
}

inline std::tuple<MathLib::Point3d, MathLib::Point3d, MathLib::Point3d>
getEdgeMiddlePoints(GeoLib::Triangle const& tri)
{
    return std::make_tuple(
        getEdgeMiddlePoint(*tri.getPoint(0), *tri.getPoint(1)),
        getEdgeMiddlePoint(*tri.getPoint(1), *tri.getPoint(2)),
        getEdgeMiddlePoint(*tri.getPoint(2), *tri.getPoint(0)));
}

/// Computes rotation matrix according to the z,x',z'' convention
inline MathLib::DenseMatrix<double, std::size_t>
getRotMat(double alpha, double beta, double gamma)
{
    MathLib::DenseMatrix<double, std::size_t> rot_mat(3,3);
    rot_mat(0,0) = cos(alpha)*cos(gamma) - sin(alpha)*cos(beta)*sin(gamma);
    rot_mat(0,1) = sin(alpha)*cos(gamma) + cos(alpha)*cos(beta)*sin(gamma);
    rot_mat(0,2) = sin(beta)*sin(gamma);
    rot_mat(1,0) = -cos(alpha)*sin(gamma) - sin(alpha)*cos(beta)*cos(gamma);
    rot_mat(1,1) = -sin(alpha)*sin(gamma) + cos(alpha)*cos(beta)*cos(gamma);
    rot_mat(1,2) = sin(beta)*cos(gamma);
    rot_mat(2,0) = sin(alpha)*sin(beta);
    rot_mat(2,1) = -cos(alpha)*sin(beta);
    rot_mat(2,2) = cos(beta);
    return rot_mat;
}

TEST(GeoLib, SurfaceIsPointInSurface)
{
    std::vector<std::function<double(double, double)>> surface_functions;
    surface_functions.emplace_back(constant);
    surface_functions.emplace_back(coscos);

    for (const auto& f : surface_functions) {
        std::random_device rd;

        std::string name("Surface");
        // generate ll and ur in random way
        std::mt19937 random_engine_mt19937(rd());
        std::normal_distribution<> normal_dist_ll(-10, 2);
        std::normal_distribution<> normal_dist_ur(10, 2);
        MathLib::Point3d ll(std::array<double,3>({{
            normal_dist_ll(random_engine_mt19937),
            normal_dist_ll(random_engine_mt19937),
            0.0}}));
        MathLib::Point3d ur(std::array<double,3>({{
            normal_dist_ur(random_engine_mt19937),
            normal_dist_ur(random_engine_mt19937),
            0.0}}));
        for (std::size_t k(0); k<3; ++k)
            if (ll[k] > ur[k])
                std::swap(ll[k], ur[k]);

        // random discretization of the domain
        std::default_random_engine re(rd());
        std::uniform_int_distribution<std::size_t> uniform_dist(2, 25);
        std::array<std::size_t,2> n_steps = {{uniform_dist(re),uniform_dist(re)}};

        std::unique_ptr<MeshLib::Mesh> sfc_mesh(
            MeshLib::MeshGenerator::createSurfaceMesh(
                name, ll, ur, n_steps, f
            )
        );

        // random rotation angles
        std::normal_distribution<> normal_dist_angles(
            0, boost::math::double_constants::two_pi);
        std::array<double,3> euler_angles = {{
            normal_dist_angles(random_engine_mt19937),
            normal_dist_angles(random_engine_mt19937),
            normal_dist_angles(random_engine_mt19937)
            }};

        MathLib::DenseMatrix<double, std::size_t> rot_mat(getRotMat(
            euler_angles[0], euler_angles[1], euler_angles[2]));

        std::vector<MeshLib::Node*> const& nodes(sfc_mesh->getNodes());
        GeoLib::rotatePoints<MeshLib::Node>(rot_mat, nodes);

        MathLib::Vector3 const normal(0,0,1.0);
        MathLib::Vector3 const surface_normal(rot_mat * normal);
        double const scaling(1e-6);
        MathLib::Vector3 const displacement(scaling * surface_normal);

        GeoLib::GEOObjects geometries;
        MeshLib::convertMeshToGeo(*sfc_mesh, geometries);

        std::vector<GeoLib::Surface*> const& sfcs(*geometries.getSurfaceVec(name));
        GeoLib::Surface const*const sfc(sfcs.front());
        std::vector<GeoLib::Point*> const& pnts(*geometries.getPointVec(name));

        double const eps(std::numeric_limits<double>::epsilon());

        // test triangle edge point of the surface triangles
        for (auto const p : pnts) {
            EXPECT_TRUE(sfc->isPntInSfc(*p, eps));
            MathLib::Point3d q(*p);
            for (std::size_t k(0); k<3; ++k)
                q[k] += displacement[k];
            EXPECT_FALSE(sfc->isPntInSfc(q, eps));
        }
        // test edge middle points of the triangles
        for (std::size_t k(0); k<sfc->getNumberOfTriangles(); ++k) {
            MathLib::Point3d p, q, r;
            std::tie(p,q,r) = getEdgeMiddlePoints(*(*sfc)[k]);
            EXPECT_TRUE(sfc->isPntInSfc(p, eps));
            EXPECT_TRUE(sfc->isPntInSfc(q, eps));
            EXPECT_TRUE(sfc->isPntInSfc(r, eps));
        }
    }
}
