/**
 * \file
 * \author Lars Bilke
 * \date   2009-10-27
 * \brief  Definition of the Vector3 class.
 *         From: http://www.strout.net/info/coding/classlib/intro.html
 *         with modifications to derive from TemplatePoint
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTOR3_H
#define VECTOR3_H

// MathLib
#include "TemplatePoint.h"

// MathLib
#include "MathTools.h"

#include <cmath>
#include <iostream>

namespace MathLib
{
/**
 * The Vector3 class defines a three-dimensional vector, with appropriate
 *	operators.  (* is cross product.)
 */
template <class T>
class TemplateVector3 : public MathLib::TemplatePoint<T>
{
public:
	TemplateVector3() : MathLib::TemplatePoint<T>() {}
	TemplateVector3(T x1, T x2, T x3) : MathLib::TemplatePoint<T>(x1, x2, x3) {}
	TemplateVector3(const MathLib::TemplatePoint<T> & rhs) :
		MathLib::TemplatePoint<T>(rhs[0], rhs[1], rhs[2])
	{}
	/** constructs a vector from the gien points */
	TemplateVector3(const MathLib::TemplatePoint<T> &a, const MathLib::TemplatePoint<T> &b) :
		MathLib::TemplatePoint<T>(b[0] - a[0], b[1] - a[1], b[2] - a[2])
	{}
	~TemplateVector3() {}

	// vector arithmetic

	TemplateVector3 operator+(const TemplateVector3 & pV) const
	{
		return TemplateVector3(this->x[0] + pV[0], this->x[1] + pV[1], this->x[2] + pV[2] );
	}

	TemplateVector3 operator-(const TemplateVector3 & pV) const
	{
		TemplateVector3 out( this->x[0] - pV[0], this->x[1] - pV[1], this->x[2] - pV[2] );
		return out;
	}

	TemplateVector3 operator-() const
	{ return TemplateVector3 (-this->x[0], -this->x[1], -this->x[2]); }

	TemplateVector3& operator+=(const TemplateVector3 & pV)
	{
		for (std::size_t i(0); i < 3; i++) this->x[i] += pV[i];
		return *this;
	}

	TemplateVector3& operator+=(const MathLib::TemplatePoint<T>& pnt)
	{
		for (std::size_t i(0); i < 3; i++) this->_x[i] += pnt[i];
		return *this;
	}

	TemplateVector3& operator-=(const TemplateVector3 & pV)
	{
		for (std::size_t i(0); i < 3; i++) this->_x[i] -= pV[i];
		return *this;
	}

	// Accessors
	T X() const { return this->_x[0]; }
	T Y() const { return this->_x[1]; }
	T Z() const { return this->_x[2]; }
	void setX(T value) { this->_x[0] = value; }
	void setY(T value) { this->_x[1] = value; }
	void setZ(T value) { this->_x[2] = value; }

	/// Dot product with another vector
	double Dot(const TemplateVector3 & pV) const
	{
		return this->_x[0] * pV[0] + this->_x[1] * pV[1] + this->_x[2] * pV[2];
	}

	/// Cross product as the multiplication operator
	TemplateVector3 operator*(const TemplateVector3 & pV) const
	{
		return TemplateVector3 (
		               this->_x[1] * pV[2] - this->_x[2] * pV[1],
		               this->_x[2] * pV[0] - this->_x[0] * pV[2],
		               this->_x[0] * pV[1] - this->_x[1] * pV[0] );
	}

	/// Cross product with another vector
	TemplateVector3 Cross( const TemplateVector3 & pV ) const
	{ return *this * pV; }

	friend double Dot( const TemplateVector3 & p1, const TemplateVector3 & p2 )
	{ return p1.Dot(p2); }

	friend TemplateVector3 Cross( const TemplateVector3 & p1, const TemplateVector3 & p2 )
	{ return p1 * p2; }

	TemplateVector3 operator*(double pR) const   // * a scalar
	{
		return TemplateVector3(this->x[0] * pR, this->x[1] * pR, this->x[2] * pR);
	}

	friend TemplateVector3 operator*(double pR, const TemplateVector3 & pV)
	{
		return TemplateVector3 ( pV[0] * pR, pV[1] * pR, pV[2] * pR );
	}

	TemplateVector3& operator*=(double pR)
	{
		for (std::size_t i(0); i < 3; i++) this->_x[i] *= pR;
		return *this;
	}

	/// Returns the squared length
	double LenSqr(void) const
	{
		return scpr<double,3> (this->getCoords (), this->getCoords ());
	}

	/// Returns the length
	double Length(void) const
	{ return sqrt(LenSqr()); }

	/// Projection (component of *this parallel to pV).
	/// Note: component perpendicular to pV is:  *this - Proj(pV)
	TemplateVector3 Proj(const TemplateVector3 & pV)
	{ TemplateVector3 out( pV * (this->Dot(pV) / pV.LenSqr()) ); return out; }

	/// Cosine of the angle between two vectors:
	double CosAng(const TemplateVector3& pV)
	{ return this->Dot(pV) / (Length() * pV.Length()); }

	/// Comparison if equal
	bool operator==( const TemplateVector3 & pV) const
	{
		return std::fabs(this->_x[0] - pV[0]) < sqrt(std::numeric_limits<double>::min()) &&
		       std::fabs(this->_x[1] - pV[1]) < sqrt(std::numeric_limits<double>::min()) &&
		       std::fabs(this->_x[2] - pV[2]) < sqrt(std::numeric_limits<double>::min());
	}

	/// Comparison if not equal
	bool operator!=( const TemplateVector3 & pV) const
	{
		return !(pV == this);
//		this->_x[0]!=pV[0] || this->_x[1]!=pV[1] || this->_x[2]!=pV[2];
	}
};

typedef TemplateVector3<double> Vector3;
}

#endif // VECTOR3_H
