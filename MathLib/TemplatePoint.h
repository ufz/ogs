/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-28
 * \brief  Definition of the TemplatePoint class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEPOINT_H_
#define TEMPLATEPOINT_H_

// STL
#include <array>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <istream>
#include <cassert>

namespace MathLib
{
/**
 * \ingroup GeoLib
 *
 * \brief class-template for points can be instantiated by a numeric type.
 * \tparam T the coordinate type
 */
template <typename T, std::size_t DIM = 3> class TemplatePoint
{
public:
    typedef T FP_T;

    /** default constructor with zero coordinates */
    TemplatePoint();

    /** constructor - constructs a TemplatePoint object
     *
     * @param x std::array containing the coordinates of the point
     */
    explicit TemplatePoint(std::array<T,DIM> const& x);

    /** virtual destructor */
    virtual ~TemplatePoint() = default;

    TemplatePoint(TemplatePoint const&) = default;
    TemplatePoint& operator=(TemplatePoint const&) = default;

    /** \brief const access operator
     *  The access to the point coordinates is like the access to a field. Code example:
     * \code
     * Point<double> point (1.0, 2.0, 3.0);
     * double sqrNrm2 = point[0] * point[0] + point[1] * point[1] + point[2] + point[2];
     * \endcode
     */
    const T& operator[] (std::size_t idx) const
    {
        assert (idx < DIM);
        return _x[idx];
    }
    /** \brief access operator (see book Effektiv C++ programmieren - subsection 1.3.2 ).
     * \sa const T& operator[] (std::size_t idx) const
     */
    T& operator[] (std::size_t idx)
    {
        return const_cast<T&> (static_cast<const TemplatePoint&> (*this)[idx]);
    }

    /** returns an array containing the coordinates of the point */
    const T* getCoords () const
    {
        return _x.data();
    }

    /** write point coordinates into stream (used from operator<<)
     * \param os a standard output stream
     */
    virtual void write (std::ostream &os) const
    {
        std::copy(_x.cbegin(), _x.cend(), std::ostream_iterator<T>(os, " "));
    }

    /** read point coordinates into stream (used from operator>>) */
    virtual void read (std::istream &is)
    {
        std::copy(std::istream_iterator<T>(is), std::istream_iterator<T>(), _x.begin());
    }

protected:
    std::array<T, DIM> _x;
};

template <typename T, std::size_t DIM>
TemplatePoint<T,DIM>::TemplatePoint() :
    _x({{0}})
{}

template <typename T, std::size_t DIM>
TemplatePoint<T,DIM>::TemplatePoint(std::array<T,DIM> const& x) :
    _x(x)
{}

/** Equality of TemplatePoint's up to an epsilon.
 */
template <typename T, std::size_t DIM>
bool operator==(TemplatePoint<T,DIM> const& a, TemplatePoint<T,DIM> const& b)
{
    T const sqr_dist(sqrDist(a,b));
    auto const eps = std::numeric_limits<T>::epsilon();
    return (sqr_dist < eps*eps);
}

template <typename T, std::size_t DIM>
bool operator< (TemplatePoint<T,DIM> const& a, TemplatePoint<T,DIM> const& b)
{
    for (std::size_t i = 0; i < DIM; ++i)
    {
        if (a[i] > b[i]) {
            return false;
        } else {
            if (a[i] < b[i]) {
                return true;
            }
        }
        // continue with next dimension, because a[0] == b[0]
    }

    // The values in all dimenisions are equal.
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
 * \f$  |a_i - b_i| > \epsilon \f$ for all coordinates \f$ 0 \le i < \textrm{DIM} \f$.
 */
template <typename T, std::size_t DIM>
bool lessEq(TemplatePoint<T, DIM> const& a, TemplatePoint<T, DIM> const& b,
        double eps = std::numeric_limits<double>::epsilon())
{
    auto coordinateIsLargerEps = [&eps](T const u, T const v) -> bool
    {
        return std::fabs(u - v) > eps * std::min(std::fabs(v), std::fabs(u)) &&
               std::fabs(u - v) > eps;
    };

    for (std::size_t i = 0; i < DIM; ++i)
    {
        // test a relative and an absolute criterion
        if (coordinateIsLargerEps(a[i], b[i]))
        {
            if (a[i] <= b[i])
                return true;
            else
                return false;
        }
        // a[i] ~= b[i] up to an epsilon. Compare next dimension.
    }

    // all coordinates are equal up to an epsilon.
    return true;
}

/** Distance between points p0 and p1 in the maximum norm. */
template <typename T>
T maxNormDist(const MathLib::TemplatePoint<T>* p0, const MathLib::TemplatePoint<T>* p1)
{
    const T x = fabs((*p1)[0] - (*p0)[0]);
    const T y = fabs((*p1)[1] - (*p0)[1]);
    const T z = fabs((*p1)[2] - (*p0)[2]);

    return std::max(x, std::max(y, z));
}

/** overload the output operator for class Point */
template <typename T, std::size_t DIM>
std::ostream& operator<< (std::ostream &os, const TemplatePoint<T,DIM> &p)
{
    p.write (os);
    return os;
}

/** overload the input operator for class Point */
template <typename T, std::size_t DIM>
std::istream& operator>> (std::istream &is, TemplatePoint<T,DIM> &p)
{
    p.read (is);
    return is;
}
} // end namespace MathLib

#endif /* TEMPLATEPOINT_H_ */
