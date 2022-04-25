/**
 * \file
 * \date   2015-01-16
 * \brief  Definition of the Point3d class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <limits>
#include <Eigen/Dense>

#include "mathlib_export.h"

namespace MathLib
{
/**
 * \ingroup MathLib
 *
 */
class Point3d
{
public:
    /** default constructor with zero coordinates */
    Point3d();

    /** constructor - constructs a Point3d object
     *
     * @param x std::array containing the coordinates of the point
     */
    explicit Point3d(std::array<double, 3> x);

    /** virtual destructor */
    virtual ~Point3d() = default;

    Point3d(Point3d const&) = default;
    Point3d& operator=(Point3d const&) = default;

    /** \brief const access operator
     *  The access to the point coordinates is like the access to a field. Code
     * example: \code Point<double> point (1.0, 2.0, 3.0); double sqrNrm2 =
     * point[0] * point[0] + point[1] * point[1] + point[2] + point[2]; \endcode
     */
    const double& operator[](std::size_t idx) const
    {
        assert(idx < 3);
        return x_[idx];
    }
    /** \brief access operator (see book Effektiv C++ programmieren -
     * subsection 1.3.2 ). \sa const T& operator[] (std::size_t idx) const
     */
    double& operator[](std::size_t idx)
    {
        return const_cast<double&>(static_cast<const Point3d&>(*this)[idx]);
    }

    /** returns an array containing the coordinates of the point */
    const double* data() const { return x_.data(); }

    double* data() { return x_.data(); }

    Eigen::Vector3d const& asEigenVector3d() const { return x_; }

private:
    Eigen::Vector3d x_;
};

inline bool operator<(Point3d const& a, Point3d const& b)
{
    return std::lexicographical_compare(a.data(), a.data() + 3, b.data(),
                                        b.data() + 3);
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
bool inline lessEq(Point3d const& a, Point3d const& b,
                   double eps = std::numeric_limits<double>::epsilon())
{
    auto absAndRelDiffLargerThanEps = [eps](double const u,
                                            double const v) -> bool
    {
        return std::abs(u - v) > eps * std::min(std::abs(v), std::abs(u)) &&
               std::abs(u - v) > eps;
    };

    return std::lexicographical_compare(
        a.data(), a.data() + 3, b.data(), b.data() + 3,
        [&absAndRelDiffLargerThanEps](auto const u, auto const v)
        {
            if (absAndRelDiffLargerThanEps(u, v))
            {
                return u <= v;
            }
            return true;
        });
}

/** overload the output operator for class Point */
inline std::ostream& operator<<(std::ostream& os, const Point3d& p)
{
    os << p[0] << " " << p[1] << " " << p[2];
    return os;
}

extern MATHLIB_EXPORT const Point3d ORIGIN;
/**
 * rotation of points
 * @param mat a rotation matrix
 * @param p   a point to be transformed
 * @return a rotated point
 */
template <typename MATRIX>
inline MathLib::Point3d operator*(MATRIX const& mat, MathLib::Point3d const& p)
{
    MathLib::Point3d new_p;
    for (std::size_t i(0); i < 3; ++i)
    {
        for (std::size_t j(0); j < 3; ++j)
        {
            new_p[i] += mat(i, j) * p[j];
        }
    }
    return new_p;
}

/** Computes the squared dist between the two points p0 and p1.
 */
double sqrDist(MathLib::Point3d const& p0, MathLib::Point3d const& p1);

/** Equality of Point3d's up to an epsilon.
 */
inline bool operator==(Point3d const& a, Point3d const& b)
{
    auto const sqr_dist(sqrDist(a, b));
    auto const eps = std::numeric_limits<double>::epsilon();
    return (sqr_dist < eps * eps);
}

/// Computes the squared distance between the orthogonal projection of the two
/// points \c p0 and \c p1 onto the \f$xy\f$-plane.
inline double sqrDist2d(MathLib::Point3d const& p0, MathLib::Point3d const& p1)
{
    return (p0[0] - p1[0]) * (p0[0] - p1[0]) +
           (p0[1] - p1[1]) * (p0[1] - p1[1]);
}

}  // end namespace MathLib
