/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <memory>
#include <numeric>
#include <random>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshQuality/ElementQualityInterface.h"
#include "MeshLib/Node.h"

class PyramidElementQuality : public ::testing::Test
{
public:
    PyramidElementQuality()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, 10);

        lengths = {1.0, 1.0, 1.0};
        n_subdivisions = {distrib(gen), distrib(gen), distrib(gen)};
    }

    std::vector<double> getElementQualityVectorFromRegularPyramidMesh(
        MeshLib::MeshQualityType const type)
    {
        std::vector<std::unique_ptr<BaseLib::ISubdivision>> vec_div;
        for (int i = 0; i < 3; ++i)
        {
            vec_div.emplace_back(
                new BaseLib::UniformSubdivision(lengths[i], n_subdivisions[i]));
        }

        std::unique_ptr<MeshLib::Mesh> mesh(
            MeshLib::MeshGenerator::generateRegularPyramidMesh(
                *vec_div[0], *vec_div[1], *vec_div[2]));
        MeshLib::ElementQualityInterface element_quality(*mesh, type);
        return element_quality.getQualityVector();
    }

    std::array<double, 3> lengths;
    std::array<int, 3> n_subdivisions;
};

TEST_F(PyramidElementQuality, ElementSize)
{
    auto const type = MeshLib::MeshQualityType::ELEMENTSIZE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularPyramidMesh(type);
    auto const expected_value =
        (1.0 / 6) / std::accumulate(n_subdivisions.begin(),
                                    n_subdivisions.end(), 1,
                                    std::multiplies<int>());
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    std::numeric_limits<double>::epsilon());
    }
}

TEST_F(PyramidElementQuality, SizeDifference)
{
    auto const type = MeshLib::MeshQualityType::SIZEDIFFERENCE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularPyramidMesh(type);
    // all elements have the same size, the quality value has to be 1.0
    auto constexpr expected_value = 1.0;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    12 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(PyramidElementQuality, EdgeRatio)
{
    auto const type = MeshLib::MeshQualityType::EDGERATIO;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularPyramidMesh(type);

    std::array const dx = {1.0 / n_subdivisions[0], 1.0 / n_subdivisions[1],
                           1.0 / n_subdivisions[2]};
    std::array const edge_lengths_top_bottom = {
        dx[0], dx[1],
        std::sqrt(std::pow(dx[0] / 2, 2) + std::pow(dx[1] / 2, 2) +
                  std::pow(dx[2] / 2, 2))};
    auto const& min_max_top_bottom = std::minmax_element(
        edge_lengths_top_bottom.begin(), edge_lengths_top_bottom.end());
    auto const expected_value_top_bottom =
        *min_max_top_bottom.first / *min_max_top_bottom.second;

    std::array const edge_lengths_left_right = {
        dx[1], dx[2],
        std::sqrt(std::pow(dx[0] / 2, 2) + std::pow(dx[1] / 2, 2) +
                  std::pow(dx[2] / 2, 2))};
    auto const& min_max_left_right = std::minmax_element(
        edge_lengths_left_right.begin(), edge_lengths_left_right.end());
    auto const expected_value_left_right =
        *min_max_left_right.first / *min_max_left_right.second;

    std::array const edge_lengths_front_back = {
        dx[0], dx[2],
        std::sqrt(std::pow(dx[0] / 2, 2) + std::pow(dx[1] / 2, 2) +
                  std::pow(dx[2] / 2, 2))};
    auto const& min_max_front_back = std::minmax_element(
        edge_lengths_front_back.begin(), edge_lengths_front_back.end());
    auto const expected_value_front_back =
        *min_max_front_back.first / *min_max_front_back.second;
    for (size_t i = 0; i < element_quality_vector.size(); i = i + 6)
    {
        ASSERT_NEAR(expected_value_top_bottom, element_quality_vector[i],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_top_bottom, element_quality_vector[i + 1],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_left_right, element_quality_vector[i + 2],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_left_right, element_quality_vector[i + 3],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_front_back, element_quality_vector[i + 4],
                    10 * std::numeric_limits<double>::epsilon());
        ASSERT_NEAR(expected_value_front_back, element_quality_vector[i + 5],
                    10 * std::numeric_limits<double>::epsilon());
    }
}
