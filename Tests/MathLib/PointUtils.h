/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <gtest/gtest.h>

#include <vector>

#include "MathLib/Point3d.h"
#include "MathLib/WeightedPoint.h"
#include "Tests/Utils.h"

namespace PointUtils
{

// Trims an array to a definite size.
template <std::size_t target_size, std::size_t source_size>
std::array<double, target_size> trimArray(
    std::array<double, source_size> const& arr)
{
    static_assert(target_size <= source_size);
    std::array<double, target_size> res{};

    std::copy(begin(arr), begin(arr) + target_size, begin(res));

    return res;
}

template <unsigned dim>
MathLib::WeightedPoint getWeightedPointOfDim(
    std::array<double, 3> const& coords)
{
    auto const coords_right_dim = trimArray<dim>(coords);

    return {coords_right_dim,
            // Weight does not matter for this unit test suite
            std::numeric_limits<double>::quiet_NaN()};
}

template <typename BulkElementType>
MathLib::WeightedPoint getWeightedPointOnFace(
    std::array<double, 3> const& face_node_natural_coords)
{
    constexpr unsigned face_dimension = BulkElementType::dimension - 1;

    return getWeightedPointOfDim<face_dimension>(face_node_natural_coords);
}

// The elements of CoordsContainer are std::array<double, 3> of coordinates.
template <unsigned dim, typename CoordsContainer>
std::vector<MathLib::WeightedPoint> toWeightedPointsOfDim(
    CoordsContainer const& coordss)
{
    std::vector<MathLib::WeightedPoint> wps;
    wps.reserve(coordss.size());

    for (auto const& coords : coordss)
    {
        wps.emplace_back(getWeightedPointOfDim<dim>(coords));
    }

    return wps;
}

// The elements of CoordsContainer are std::array<double, 3> of coordinates.
template <typename BulkElementType, typename CoordsContainer>
std::vector<MathLib::WeightedPoint> toWeightedPointsOnFace(
    CoordsContainer const& coordss)
{
    return toWeightedPointsOfDim<BulkElementType::dimension - 1>(coordss);
}

inline std::array<double, 3> getCoords(MathLib::Point3d const& p)
{
    return {p[0], p[1], p[2]};
}

template <typename T>
struct PrintableRef
{
    PrintableRef(T const& data) : data_{data} {}

    friend std::ostream& operator<<(std::ostream& os,
                                    PrintableRef<T> const& ref)
    {
        return os << ref.data_;
    }

private:
    T const& data_;
};

template <typename T, std::size_t N>
struct PrintableRef<std::array<T, N>>
{
private:
    using Data = std::array<T, N>;

public:
    PrintableRef(Data const& data) : data_{data} {}

    friend std::ostream& operator<<(std::ostream& os,
                                    PrintableRef<Data> const& ref)
    {
        auto const& arr = ref.data_;
        os << '[';
        for (std::size_t i = 0; i < N; ++i)
        {
            if (i != 0)
            {
                os << ", ";
            }
            os << arr[i];
        }
        os << ']';
        return os;
    }

private:
    Data const& data_;
};

#ifdef OGS_HAVE_CONCEPTS
template <typename T>
concept PointLike =
    std::same_as<T, MathLib::Point3d> || std::same_as<T, std::array<double, 3>>;
#endif

struct IsNear
{
    template <OGS_USE_CONCEPT(PointLike) PointA,
              OGS_USE_CONCEPT(PointLike) PointB>
    testing::AssertionResult operator()(const char* a_expr, const char* b_expr,
                                        const char* /*abstol_expr*/,
                                        PointA const& a, PointB const& b,
                                        double const abstol) const
    {
        for (std::size_t d : {0, 1, 2})
        {
            auto const diff = std::abs(a[d] - b[d]);
            if (std::isnan(diff) || diff > abstol)
            {
                return testing::AssertionFailure()
                       << a_expr << " and " << b_expr << " differ by " << diff
                       << " > " << abstol << " in coordinate #" << d << '\n'
                       << a_expr << " evaluates to " << PrintableRef{a} << '\n'
                       << b_expr << " evaluates to " << PrintableRef{b};
            }
        }

        return testing::AssertionSuccess();
    }

    template <OGS_USE_CONCEPT(PointLike) PointA,
              OGS_USE_CONCEPT(PointLike) PointB>
    testing::AssertionResult operator()(const char* a_expr, const char* b_expr,
                                        PointA const& a, PointB const& b) const
    {
        return (*this)(a_expr, b_expr, "NO EXPR AVAIL", a, b,
                       std::numeric_limits<double>::epsilon());
    }
};
}
