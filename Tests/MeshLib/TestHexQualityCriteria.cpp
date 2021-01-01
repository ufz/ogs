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

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshQuality/ElementQualityInterface.h"

class HexElementQuality : public ::testing::Test
{
public:
    HexElementQuality()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(1, 10);

        lengths = {1, 1, 1};
        n_subdivisions = {distrib(gen), distrib(gen), distrib(gen)};
    }

    std::vector<double> getElementQualityVectorFromRegularHexMesh(
        MeshLib::MeshQualityType const type)
    {
        std::vector<std::unique_ptr<BaseLib::ISubdivision>> vec_div;
        for (int i = 0; i < 3; ++i)
        {
            vec_div.emplace_back(
                new BaseLib::UniformSubdivision(lengths[i], n_subdivisions[i]));
        }

        std::unique_ptr<MeshLib::Mesh> mesh(
            MeshLib::MeshGenerator::generateRegularHexMesh(
                *vec_div[0], *vec_div[1], *vec_div[2]));
        MeshLib::ElementQualityInterface element_quality(*mesh, type);
        return element_quality.getQualityVector();
    }

    std::array<int, 3> lengths;
    std::array<int, 3> n_subdivisions;
};

TEST_F(HexElementQuality, ElementSize)
{
    auto const type = MeshLib::MeshQualityType::ELEMENTSIZE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularHexMesh(type);
    auto const expected_value =
        1.0 / std::accumulate(n_subdivisions.begin(), n_subdivisions.end(), 1,
                              std::multiplies<int>());
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    std::numeric_limits<double>::epsilon());
    }
}

TEST_F(HexElementQuality, SizeDifference)
{
    auto const type = MeshLib::MeshQualityType::SIZEDIFFERENCE;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularHexMesh(type);
    // all elements have the same size, the quality value has to be 1.0
    auto constexpr expected_value = 1.0;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(HexElementQuality, EdgeRatio)
{
    auto const type = MeshLib::MeshQualityType::EDGERATIO;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularHexMesh(type);
    auto const& min_max =
        std::minmax_element(n_subdivisions.begin(), n_subdivisions.end());
    auto const expected_value = double(*min_max.first) / *min_max.second;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    10 * std::numeric_limits<double>::epsilon());
    }
}

TEST_F(HexElementQuality, EquiAngleSkew)
{
    auto const type = MeshLib::MeshQualityType::EQUIANGLESKEW;
    auto const element_quality_vector =
        getElementQualityVectorFromRegularHexMesh(type);
    auto const expected_value = 0.0;
    for (auto const element_quality : element_quality_vector)
    {
        ASSERT_NEAR(expected_value, element_quality,
                    2*std::numeric_limits<double>::epsilon());
    }
}
