/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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

class TriElementQuality : public ::testing::Test
{
public:
    TriElementQuality()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, 10);

        lengths = {1, 1};
        n_subdivisions = {distrib(gen), distrib(gen)};
    }

    std::vector<double> getElementQualityVectorFromRegularTriMesh(
        MeshLib::MeshQualityType const type)
    {
        std::vector<std::unique_ptr<BaseLib::ISubdivision>> vec_div;
        for (int i = 0; i < 2; ++i)
        {
            vec_div.emplace_back(
                new BaseLib::UniformSubdivision(lengths[i], n_subdivisions[i]));
        }

        std::unique_ptr<MeshLib::Mesh> mesh(
            MeshLib::MeshGenerator::generateRegularTriMesh(*vec_div[0],
                                                           *vec_div[1]));
        MeshLib::ElementQualityInterface element_quality(*mesh, type);
        return element_quality.getQualityVector();
    }

    std::array<int, 2> lengths;
    std::array<int, 2> n_subdivisions;
};

TEST_F(TriElementQuality, ElementSize)
{
    auto const type = MeshLib::MeshQualityType::ELEMENTSIZE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTriMesh(type);
    auto const expected_value =
        0.5 / std::accumulate(n_subdivisions.begin(), n_subdivisions.end(), 1,
                              std::multiplies<int>());
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    std::numeric_limits<double>::epsilon());
    }
}

TEST_F(TriElementQuality, SizeDifference)
{
    auto const type = MeshLib::MeshQualityType::SIZEDIFFERENCE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTriMesh(type);
    // all elements have the same size, the quality value has to be 1.0
    auto constexpr expected_value = 1.0;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(TriElementQuality, EdgeRatio)
{
    auto const type = MeshLib::MeshQualityType::EDGERATIO;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTriMesh(type);
    auto const& min_max =
        std::minmax_element(n_subdivisions.begin(), n_subdivisions.end());
    auto const expected_value =
        double(*min_max.first) / std::sqrt(*min_max.first * *min_max.first +
                                           *min_max.second * *min_max.second);
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(TriElementQuality, EquiAngleSkew)
{
    using namespace boost::math::double_constants;
    auto const type = MeshLib::MeshQualityType::EQUIANGLESKEW;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularTriMesh(type);
    // all triangles are right-angled triangles
    auto const hypothenuse = std::sqrt(std::pow(1.0 / n_subdivisions[0], 2) +
                                       std::pow(1.0 / n_subdivisions[1], 2));
    std::array const angles = {
        std::asin((1.0 / n_subdivisions[0]) / hypothenuse),
        std::asin((1.0 / n_subdivisions[1]) / hypothenuse), half_pi};
    auto const& min_max = std::minmax_element(angles.begin(), angles.end());
    auto const expected_value =
        std::max((*min_max.second - third_pi) / two_thirds_pi,
                 (third_pi - *min_max.first) / third_pi);

    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}
