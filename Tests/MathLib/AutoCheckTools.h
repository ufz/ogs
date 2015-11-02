/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef TESTS_MATHLIB_AUTOCHECKTOOLS_H_
#define TESTS_MATHLIB_AUTOCHECKTOOLS_H_

#include "autocheck/autocheck.hpp"

#include "MathLib/Point3d.h"

namespace autocheck
{
template <typename T, std::size_t N, typename Gen = generator<T>>
struct randomTupleGenerator
{
    Gen source;

    using result_type = std::array<T, N>;

    result_type operator()(std::size_t size = 0)
    {
        result_type rv;
        std::generate(rv.begin(), rv.end(), fix(size, source));
        return rv;
    }
};

enum class CartesianPlane
{
	XY = 2,  // The values are used for indexing.
	YZ = 0,
	ZX = 1
};

enum class CartesianAxes
{
	X = 0,  // The values are used for indexing.
	Y = 1,
	Z = 2
};

template <enum CartesianPlane P, typename T,
          typename Gen = randomTupleGenerator<T, 2>>
struct tripleInPlaneGenerator
{
	Gen source;

	using result_type = std::array<T, 3>;

	result_type operator()(std::size_t size = 0)
	{
		typename Gen::result_type const tuple = fix(size, source)();
		result_type rv;
		rv.fill(T());	// fill with default value for type T.
		int j = 0;	// running over the tuple elements
		for (int i = 0; i < 3; ++i) // i is running over the rv elems.
		{
			if (i == static_cast<int>(P))
				continue;
			rv[i] = tuple[j++];
		}

		return rv;
    }
};

template <enum CartesianAxes A, typename T,
          typename Gen = randomTupleGenerator<T, 1>>
struct tripleOnAxisGenerator
{
	Gen source;

	using result_type = std::array<T, 3>;

	result_type operator()(std::size_t size = 0)
	{
		typename Gen::result_type const tuple = fix(size, source)();
		result_type rv;
		rv.fill(T());	// fill with default value for type T.
		int j = 0;	// running over the tuple elements
		for (int i = 0; i < 3; ++i) // i is running over the rv elems.
		{
			if (i == static_cast<int>(A))
				rv[i] = tuple[j++];
		}

		return rv;
    }
};
/// Generates non-negative integers from 0 to given maximum dimension DIM
/// independent of size.
template <typename T, T DIM, typename Gen = generator<T>>
struct randomCoordinateIndexGenerator
{
    Gen source;
    using result_type = T;

    result_type operator()(std::size_t)
    {
        return source() % DIM;
    }
};

static double absoluteValue(double&& v, std::size_t)
{
    return std::abs(v);
}

/// Generates values of T starting with given maximum value (default 1) and
/// progressively smaller as size increases.
template <typename T, typename Gen = generator<T>>
struct progressivelySmallerGenerator
{
    Gen source;
    T _max_value;

    using result_type = T;

    progressivelySmallerGenerator(T max_value = T{1}) : _max_value(max_value)
    {
    }

    result_type operator()(std::size_t size = 0)
    {
        return _max_value * fix(1, source)() / (size + 1);
    }
};

}  // namespace autocheck
#endif  // TESTS_MATHLIB_AUTOCHECKTOOLS_H_
