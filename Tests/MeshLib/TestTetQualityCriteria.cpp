/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshQuality/ElementQualityInterface.h"
#include "MeshLib/Node.h"
#include "gtest/gtest.h"

class TetElementQuality : public ::testing::Test
{
public:
    TetElementQuality()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, 10);

        lengths = {1, 1, 1};
        n_subdivisions = {distrib(gen), distrib(gen), distrib(gen)};
    }

    std::vector<double> getElementQualityVectorFromRegularTetMesh(
        MeshLib::MeshQualityType const type)
    {
        std::vector<std::unique_ptr<BaseLib::ISubdivision>> vec_div;
        for (int i = 0; i < 3; ++i)
        {
            vec_div.emplace_back(
                new BaseLib::UniformSubdivision(lengths[i], n_subdivisions[i]));
        }

        std::unique_ptr<MeshLib::Mesh> mesh(
            MeshLib::MeshGenerator::generateRegularTetMesh(
                *vec_div[0], *vec_div[1], *vec_div[2]));
        MeshLib::ElementQualityInterface element_quality(*mesh, type);
        return element_quality.getQualityVector();
    }

    std::array<int, 3> lengths;
    std::array<int, 3> n_subdivisions;
};

TEST_F(TetElementQuality, ElementSize)
{
    auto const type = MeshLib::MeshQualityType::ELEMENTSIZE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTetMesh(type);
    auto const expected_value =
        1.0 / (6 * std::accumulate(n_subdivisions.begin(), n_subdivisions.end(),
                                   1, std::multiplies<int>()));
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    std::numeric_limits<double>::epsilon());
    }
}

TEST_F(TetElementQuality, SizeDifference)
{
    auto const type = MeshLib::MeshQualityType::SIZEDIFFERENCE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTetMesh(type);
    // all elements have the same size, the quality value has to be 1.0
    auto constexpr expected_value = 1.0;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(TetElementQuality, EdgeRatio)
{
    auto const type = MeshLib::MeshQualityType::EDGERATIO;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTetMesh(type);
    std::array const d = {1.0 / n_subdivisions[0], 1.0 / n_subdivisions[1],
                          1.0 / n_subdivisions[2]};
    std::array const d2 = {std::pow(d[0], 2), std::pow(d[1], 2),
                           std::pow(d[2], 2)};
    auto const& min_max = std::minmax_element(d.begin(), d.end());
    auto const diagonal_length = std::sqrt(d2[0] + d2[1] + d2[2]);
    auto const expected_value_tet0 = *min_max.first / diagonal_length;
    auto const min_x_z = std::min(d[0], d[2]);
    auto const expected_value_tet1 = min_x_z / diagonal_length;
    auto const expected_value_tet2 =
        *min_max.first /
        std::max({std::sqrt(d2[0] + d2[1]), std::sqrt(d2[1] + d2[2]),
                  std::sqrt(d2[0] + d2[2])});
    auto const expected_value_tet3 = expected_value_tet2;
    auto const expected_value_tet4 = expected_value_tet1;
    auto const expected_value_tet5 = expected_value_tet0;

    for (int i = 0; i < element_quality_vector.size(); i = i + 6)
    {
        ASSERT_NEAR(expected_value_tet0, element_quality_vector[i],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_tet1, element_quality_vector[i + 1],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_tet2, element_quality_vector[i + 2],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_tet3, element_quality_vector[i + 3],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_tet4, element_quality_vector[i + 4],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_tet5, element_quality_vector[i + 5],
                    10 * std::numeric_limits<double>::epsilon());
    }
}

// compute angle formed by (n0, n1, n2)
double getAngleFromTriangle(Eigen::Vector3d n0, Eigen::Vector3d n1,
                            Eigen::Vector3d n2)
{
    auto const u = n1 - n0;
    auto const v = n1 - n2;
    return std::acos(((u.transpose() * v) / (u.norm() * v.norm()))(0, 0));
}

// get min and max angle from triangle
std::pair<double, double> getMinMaxAngleFromTriangle(Eigen::Vector3d n0,
                                                     Eigen::Vector3d n1,
                                                     Eigen::Vector3d n2)
{
    std::array const angles = {getAngleFromTriangle(n0, n1, n2),
                               getAngleFromTriangle(n1, n2, n0),
                               getAngleFromTriangle(n2, n0, n1)};
    auto const pair = std::minmax_element(angles.begin(), angles.end());
    return {*pair.first, *pair.second};
}

double computeCriterionForTet(Eigen::Vector3d n0, Eigen::Vector3d n1,
                              Eigen::Vector3d n2, Eigen::Vector3d n3)
{
    using namespace boost::math::double_constants;
    auto const alpha0 = getMinMaxAngleFromTriangle(n0, n2, n1);
    auto const alpha1 = getMinMaxAngleFromTriangle(n0, n1, n3);
    auto const alpha2 = getMinMaxAngleFromTriangle(n1, n2, n3);
    auto const alpha3 = getMinMaxAngleFromTriangle(n2, n0, n3);
    auto const min =
        std::min({alpha0.first, alpha1.first, alpha2.first, alpha3.first});
    auto const max =
        std::max({alpha0.second, alpha1.second, alpha2.second, alpha3.second});
    return std::max((max - third_pi) / two_thirds_pi,
                    (third_pi - min) / third_pi);
}

TEST_F(TetElementQuality, EquiAngleSkew)
{
    auto const type = MeshLib::MeshQualityType::EQUIANGLESKEW;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTetMesh(type);
    std::array d = {1.0 / n_subdivisions[0], 1.0 / n_subdivisions[1],
                    1.0 / n_subdivisions[2]};
    // 6 tets are within a hexahedron, the hex nodes are
    Eigen::Vector3d const n0(0, 0, 0);
    Eigen::Vector3d const n1(d[0], 0, 0);
    Eigen::Vector3d const n2(0, d[1], 0);
    Eigen::Vector3d const n3(d[0], d[1], 0);
    Eigen::Vector3d const n4(0, 0, d[2]);
    Eigen::Vector3d const n5(d[0], 0, d[2]);
    Eigen::Vector3d const n6(0, d[1], d[2]);
    Eigen::Vector3d const n7(d[0], d[1], d[2]);
    // tet 0 consist of nodes n0, n3, n2, n4
    auto const expected_value_tet0 = computeCriterionForTet(n0, n3, n2, n4);
    // tet 1 consist of nodes n3, n2, n4, n7
    auto const expected_value_tet1 = computeCriterionForTet(n3, n2, n4, n7);
    // tet 2 consist of nodes n2, n4, n7, n6
    auto const expected_value_tet2 = computeCriterionForTet(n2, n4, n7, n6);
    // tet 3 consist of nodes n0, n1, n3, n5
    auto const expected_value_tet3 = computeCriterionForTet(n0, n1, n3, n5);
    // tet 4 consist of nodes n0, n3, n4, n5
    auto const expected_value_tet4 = computeCriterionForTet(n0, n3, n4, n5);
    // tet 5 consist of nodes n3, n4, n5, n7
    auto const expected_value_tet5 = computeCriterionForTet(n3, n4, n5, n7);

    auto constexpr eps = 11 * std::numeric_limits<double>::epsilon();
    for (int i = 0; i < element_quality_vector.size(); i = i + 6)
    {
        ASSERT_NEAR(expected_value_tet0, element_quality_vector[i], eps);
        ASSERT_NEAR(expected_value_tet1, element_quality_vector[i + 1], eps);
        ASSERT_NEAR(expected_value_tet2, element_quality_vector[i + 2], eps);
        ASSERT_NEAR(expected_value_tet3, element_quality_vector[i + 3], eps);
        ASSERT_NEAR(expected_value_tet4, element_quality_vector[i + 4], eps);
        ASSERT_NEAR(expected_value_tet5, element_quality_vector[i + 5], eps);
    }
}
