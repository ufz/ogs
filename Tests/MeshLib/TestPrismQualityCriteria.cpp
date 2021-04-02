/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>
#include <numeric>
#include <random>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshQuality/ElementQualityInterface.h"
#include "MeshLib/Node.h"
#include "gtest/gtest.h"

class PrismElementQuality : public ::testing::Test
{
public:
    PrismElementQuality()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, 10);

        lengths = {1.0, 1.0, 1.0};
        n_subdivisions = {distrib(gen), distrib(gen), distrib(gen)};
    }

    std::vector<double> getElementQualityVectorFromRegularPrismMesh(
        MeshLib::MeshQualityType const type)
    {
        std::unique_ptr<MeshLib::Mesh> mesh(
            MeshLib::MeshGenerator::generateRegularPrismMesh(
                lengths[0], lengths[1], lengths[2], n_subdivisions[0],
                n_subdivisions[1], n_subdivisions[2]));
        MeshLib::ElementQualityInterface element_quality(*mesh, type);
        return element_quality.getQualityVector();
    }

    std::array<double, 3> lengths;
    std::array<int, 3> n_subdivisions;
};

TEST_F(PrismElementQuality, ElementSize)
{
    auto const type = MeshLib::MeshQualityType::ELEMENTSIZE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularPrismMesh(type);
    auto const expected_value =
        0.5 / std::accumulate(n_subdivisions.begin(), n_subdivisions.end(), 1,
                              std::multiplies<int>());
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    std::numeric_limits<double>::epsilon());
    }
}

TEST_F(PrismElementQuality, SizeDifference)
{
    auto const type = MeshLib::MeshQualityType::SIZEDIFFERENCE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularPrismMesh(type);
    // all elements have the same size, the quality value has to be 1.0
    auto constexpr expected_value = 1.0;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(PrismElementQuality, EdgeRatio)
{
    auto const type = MeshLib::MeshQualityType::EDGERATIO;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularPrismMesh(type);
    std::array const dx = {1.0 / n_subdivisions[0], 1.0 / n_subdivisions[1],
                           1.0 / n_subdivisions[2]};
    std::array const edge_lengths = {
        dx[0], dx[1], dx[2],
        std::sqrt(std::pow(dx[0], 2) + std::pow(dx[1], 2))};
    auto const& min_max =
        std::minmax_element(edge_lengths.begin(), edge_lengths.end());
    auto const expected_value = *min_max.first / *min_max.second;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(PrismElementQuality, EquiAngleSkew)
{
    auto const type = MeshLib::MeshQualityType::EQUIANGLESKEW;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularPrismMesh(type);
    auto const expected_value = 0.0;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    2 * std::numeric_limits<double>::epsilon());
    }
}
