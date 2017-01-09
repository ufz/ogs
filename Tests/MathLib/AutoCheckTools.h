/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <Eigen/Dense>
#include <autocheck/autocheck.hpp>

#include "MathLib/Point3d.h"

namespace autocheck
{

/// The autocheck standard generator for floating point values returns uniform
/// distributed values. Autocheck uses the size argument of the operator() to
/// determine the interval for the uniform distribution. This class template
/// fixes the size to one and transforms a generated random value (within the
/// interval [-1,1]) linearly to a value in the interval [a,b].
template <typename T, typename Gen = generator<T>>
struct IntervalGenerator
{
    /// Construtor initializing the slope and the \f$y\f$-intercept deploying
    /// lower bound \f$a\f$ and upper bound \f$b\f$ of the interval.
    IntervalGenerator(T a, T b)
        : _m((b-a)/2), _n((b+a)/2)
    {}

    // parameters for the interval mapping [-1,1] -> [a,b],
    // y = _m * x + _n
    T _m{1};
    T _n{0};

    Gen generator;

    using result_type = T;

    result_type intervalMap(T val) const
    {
        return _m * val + _n;
    }

    result_type operator()(std::size_t /*size*/ = 0)
    {
        return intervalMap(fix(1, generator)());
    }
};

template <typename T, typename Gen = IntervalGenerator<T>>
struct IntervalTupleGenerator
{
    IntervalTupleGenerator(Gen& ig0, Gen& ig1, Gen& ig2)
        : x_gen(ig0), y_gen(ig1), z_gen(ig2)
    {}

    Gen x_gen, y_gen, z_gen;

    using result_type = std::array<T, 3>;

    result_type operator()(std::size_t /*size*/ = 0)
    {
        return {{ x_gen(), y_gen(), z_gen() }};
    }
};

/// Generator for MxN fixed size eigen matrices with underlying type T.
template <typename T, int M, int N, typename Gen = generator<T>>
struct randomEigenMatrixGenerator
{
    Gen source;

    using result_type = Eigen::Matrix<T, M, N>;

    result_type operator()(std::size_t size = 0)
    {
        result_type rv;
        std::generate_n(rv.data(), M*N, fix(size, source));
        return rv;
    }
};

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

inline double absoluteValue(double&& v, std::size_t)
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
