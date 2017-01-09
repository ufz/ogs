/**
 * \file
 * \author Lars Bilke
 * \date   2009-10-27
 * \brief  Definition of the Vector3 class.
 *         From: http://www.strout.net/info/coding/classlib/intro.html
 *         with modifications to derive from TemplatePoint
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTOR3_H
#define VECTOR3_H

#include <cmath>

#include "TemplatePoint.h"
#include "MathTools.h"

namespace MathLib
{
/**
 * The Vector3 class defines a three-dimensional vector, with appropriate
 * operators.
 */
template <class T>
class TemplateVector3 : public MathLib::TemplatePoint<T>
{
public:
    /**
     * Default constructor. All coordinates are set to zero.
     */
    TemplateVector3() = default;

    TemplateVector3(T x0, T x1, T x2)
    {
        this->_x[0] = x0;
        this->_x[1] = x1;
        this->_x[2] = x2;
    }

    /**
     * Copy constructor.
     */
    TemplateVector3(TemplateVector3<T> const& /* v */) = default;
    TemplateVector3<T>& operator=(TemplateVector3<T> const& /* v */) = default;

    /**
     * Construct Vector3 from TemplatePoint.
     */
    TemplateVector3(TemplatePoint<T,3> const& p) :
        TemplatePoint<T>(p)
    {}

    /** Constructs the vector \f$v=(b-a)\f$ from the given points,
     * which starts in point \f$a\f$ and ends in point \f$b\f$
     */
    TemplateVector3(const MathLib::TemplatePoint<T> &a, const MathLib::TemplatePoint<T> &b) :
        MathLib::TemplatePoint<T>()
    {
        this->_x[0] = b[0] - a[0];
        this->_x[1] = b[1] - a[1];
        this->_x[2] = b[2] - a[2];
    }

    // vector arithmetic
    TemplateVector3 operator+(TemplateVector3 const& v) const
    {
        return TemplateVector3(this->_x[0]+v[0], this->_x[1]+v[1], this->_x[2]+v[2]);
    }

    TemplateVector3 operator-(TemplateVector3 const& v) const
    {
        return TemplateVector3(this->_x[0]-v[0], this->_x[1]-v[1], this->_x[2]-v[2]);
    }

    TemplateVector3& operator+=(TemplateVector3 const& v)
    {
        for (std::size_t i(0); i < 3; i++) this->_x[i] += v[i];
        return *this;
    }

    TemplateVector3& operator-=(const TemplateVector3 & pV)
    {
        for (std::size_t i(0); i < 3; i++) this->_x[i] -= pV[i];
        return *this;
    }

    TemplateVector3& operator*=(double s)
    {
        for (std::size_t i(0); i < 3; i++)
            this->_x[i] *= s;
        return *this;
    }

    /**
     * After applying the normalize operator to the vector its length is 1.0.
     */
    void normalize()
    {
        const double s(1/getLength());
        for (std::size_t i(0); i < 3; i++)
            this->_x[i] *= s;
    }

    /// Returns a normalized version of this vector
    TemplateVector3<double> getNormalizedVector() const
    {
        if (getSqrLength() == 0)
            return TemplateVector3<double>(0,0,0);
        TemplateVector3<double> norm_vec (this->_x[0], this->_x[1], this->_x[2]);
        norm_vec.normalize();
        return norm_vec;
    }

    /// Returns the squared length
    double getSqrLength(void) const
    {
        return this->_x[0]*this->_x[0] + this->_x[1]*this->_x[1] + this->_x[2]*this->_x[2];
    }

    /// Returns the length
    double getLength(void) const
    {
        return sqrt(getSqrLength());
    }

    /** scalarProduct, implementation of scalar product,
     * sometimes called dot or inner product.
     */
    template <typename T1>
    friend T1 scalarProduct(TemplateVector3<T1> const& v, TemplateVector3<T1> const& w);

    /** crossProduct: implementation of cross product,
     * sometimes called outer product.
     */
    template <typename T1>
    friend TemplateVector3<T1> crossProduct(
        TemplateVector3<T1> const& v,
        TemplateVector3<T1> const& w);

    /**  multiplication with a scalar s */
    template <typename T1>
    friend     TemplateVector3<T1> operator*(
        TemplateVector3<T1> const& v,
        double s);
    template <typename T1>
    friend     TemplateVector3<T1> operator*(
        double s,
        TemplateVector3<T1> const& v);
};

template <typename T>
T scalarProduct(TemplateVector3<T> const& v, TemplateVector3<T> const& w)
{
    return v._x[0] * w._x[0] + v._x[1] * w._x[1] + v._x[2] * w._x[2];
}

template <typename T1>
TemplateVector3<T1> crossProduct(
        TemplateVector3<T1> const& v,
        TemplateVector3<T1> const& w)
{
    return TemplateVector3<T1>(
            v._x[1] * w._x[2] - v._x[2] * w._x[1],
            v._x[2] * w._x[0] - v._x[0] * w._x[2],
            v._x[0] * w._x[1] - v._x[1] * w._x[0]);
}

template <typename T1> TemplateVector3<T1> operator*(
        TemplateVector3<T1> const& v,
        double s)
{
    return TemplateVector3<T1>(v[0] * s, v[1] * s, v[2] * s);
}

template <typename T1> TemplateVector3<T1> operator*(
        double s,
        TemplateVector3<T1> const& v)
{
    return v * s;
}

typedef TemplateVector3<double> Vector3;

/// Calculates the scalar triple (u x v) . w
double scalarTriple(MathLib::Vector3 const& u, MathLib::Vector3 const& v,
                    MathLib::Vector3 const& w);

}

#endif // VECTOR3_H
