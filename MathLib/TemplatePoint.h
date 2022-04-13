/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-28
 * \brief  Definition of the TemplatePoint class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "mathlib_export.h"

// STL
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <istream>
#include <iterator>
#include <utility>

namespace MathLib
{
/**
 * \ingroup GeoLib
 *
 * \brief class-template for points can be instantiated by a numeric type.
 * \tparam T the coordinate type
 */
template <typename T>
class TemplatePoint
{
public:
    /** default constructor with zero coordinates */
    TemplatePoint();

    /** constructor - constructs a TemplatePoint object
     *
     * @param x std::array containing the coordinates of the point
     */
    explicit TemplatePoint(std::array<T, 3> x);

    /** virtual destructor */
    virtual ~TemplatePoint() = default;

    TemplatePoint(TemplatePoint const&) = default;
    TemplatePoint& operator=(TemplatePoint const&) = default;

    /** \brief const access operator
     *  The access to the point coordinates is like the access to a field. Code
     * example: \code Point<double> point (1.0, 2.0, 3.0); double sqrNrm2 =
     * point[0] * point[0] + point[1] * point[1] + point[2] + point[2]; \endcode
     */
    const T& operator[](std::size_t idx) const
    {
        assert(idx < 3);
        return x_[idx];
    }
    /** \brief access operator (see book Effektiv C++ programmieren -
     * subsection 1.3.2 ). \sa const T& operator[] (std::size_t idx) const
     */
    T& operator[](std::size_t idx)
    {
        return const_cast<T&>(static_cast<const TemplatePoint&>(*this)[idx]);
    }

    /** returns an array containing the coordinates of the point */
    const T* getCoords() const { return x_.data(); }

    T* getCoords() { return x_.data(); }

private:
    std::array<T, 3> x_;
};

template <typename T>
TemplatePoint<T>::TemplatePoint() : x_({{0}})
{
}

template <typename T>
TemplatePoint<T>::TemplatePoint(std::array<T, 3> x) : x_(std::move(x))
{
}

/** Equality of TemplatePoint's up to an epsilon.
 */
template <typename T>
bool operator==(TemplatePoint<T> const& a, TemplatePoint<T> const& b)
{
    T const sqr_dist(sqrDist(a, b));
    auto const eps = std::numeric_limits<T>::epsilon();
    return (sqr_dist < eps * eps);
}

template <typename T>
bool operator<(TemplatePoint<T> const& a, TemplatePoint<T> const& b)
{
    for (std::size_t i = 0; i < 3; ++i)
    {
        if (a[i] > b[i])
        {
            return false;
        }
        if (a[i] < b[i])
        {
            return true;
        }

        // continue with next dimension, because a[0] == b[0]
    }

    // The values in all dimensions are equal.
    return false;
}

/**
 * Lexicographic comparison of points taking an epsilon into account.
 *
 * @param a first input point.
 * @param b second input point.
 * @param eps tolerance used in comparison of coordinates.
 *
 * @return true, if a is smaller then or equal to b according to the following
 * test \f$ |a_i - b_i| > \epsilon \cdot \min (|a_i|, |b_i|) \f$ \b and
 * \f$  |a_i - b_i| > \epsilon \f$ for all coordinates \f$ 0 \le i < 3 \f$.
 */
template <typename T>
bool lessEq(TemplatePoint<T> const& a, TemplatePoint<T> const& b,
            double eps = std::numeric_limits<double>::epsilon())
{
    auto coordinateIsLargerEps = [&eps](T const u, T const v) -> bool
    {
        return std::abs(u - v) > eps * std::min(std::abs(v), std::abs(u)) &&
               std::abs(u - v) > eps;
    };

    for (std::size_t i = 0; i < 3; ++i)
    {
        // test a relative and an absolute criterion
        if (coordinateIsLargerEps(a[i], b[i]))
        {
            return static_cast<bool>(a[i] <= b[i]);
        }
        // a[i] ~= b[i] up to an epsilon. Compare next dimension.
    }

    // all coordinates are equal up to an epsilon.
    return true;
}

/** overload the output operator for class Point */
template <typename T>
std::ostream& operator<<(std::ostream& os, const TemplatePoint<T>& p)
{
    os << p[0] << " " << p[1] << " " << p[2];
    return os;
}
}  // end namespace MathLib
