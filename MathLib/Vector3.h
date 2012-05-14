/**
 * \file Vector3.h
 * 27/10/2009 LB Initial implementation
 * From: http://www.strout.net/info/coding/classlib/intro.html
 * with modifications to derive from TemplatePoint
 */

#ifndef VECTOR3_H
#define VECTOR3_H

// GEO
#include "TemplatePoint.h"

// MathLib
#include "MathTools.h"

#include <iostream>
#include <cmath>

namespace MathLib {

/**
 * The Vector class defines a three-dimensional vector, with appropriate
 *	operators.  (* is cross product.)
 */
template <class T>
class TemplateVector : public GeoLib::TemplatePoint<T>
{
public:
	TemplateVector() : GeoLib::TemplatePoint<T>() {};
	TemplateVector(T x1, T x2, T x3) : GeoLib::TemplatePoint<T>(x1, x2, x3) {};
 	TemplateVector(const GeoLib::TemplatePoint<T> & rhs) :
 		GeoLib::TemplatePoint<T>(rhs[0], rhs[1], rhs[2])
 	{}
 	/** constructs a vector from the gien points */
 	TemplateVector(const GeoLib::TemplatePoint<T> &a, const GeoLib::TemplatePoint<T> &b) :
 	 	GeoLib::TemplatePoint<T>(b[0]-a[0], b[1]-a[1], b[2]-a[2])
 	{}
	~TemplateVector() {};

	// vector arithmetic

	TemplateVector operator+(const TemplateVector & pV) const
	{
		return TemplateVector(this->x[0]+pV[0], this->x[1]+pV[1], this->x[2]+pV[2] );
	}

	TemplateVector operator-(const TemplateVector & pV) const {
		TemplateVector out( this->x[0]-pV[0], this->x[1]-pV[1], this->x[2]-pV[2] );
		return out;
	}

	TemplateVector operator-() const
	{ return TemplateVector (-this->x[0], -this->x[1], -this->x[2]); }

	TemplateVector& operator+=(const TemplateVector & pV) {
		for (size_t i(0); i<3; i++) this->x[i]+=pV[i];
		return *this;
	}

	TemplateVector& operator+=(const GeoLib::TemplatePoint<T>& pnt)
	{
		for (size_t i(0); i<3; i++) this->_x[i] += pnt[i];
		return *this;
	}

	TemplateVector& operator-=(const TemplateVector & pV)
	{
		for (size_t i(0); i<3; i++) this->_x[i] -= pV[i];
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
	double Dot(const TemplateVector & pV) const
	{
		return this->_x[0]*pV[0] + this->_x[1]*pV[1] + this->_x[2]*pV[2];
	}

	/// Cross product as the multiplication operator
	TemplateVector operator*(const TemplateVector & pV) const {
		return TemplateVector (
				this->_x[1]*pV[2]-this->_x[2]*pV[1],
				this->_x[2]*pV[0]-this->_x[0]*pV[2],
				this->_x[0]*pV[1]-this->_x[1]*pV[0] );
	}

	/// Cross product with another vector
	TemplateVector Cross( const TemplateVector & pV ) const
	{ return *this * pV; }

	friend double Dot( const TemplateVector & p1, const TemplateVector & p2 )
	{ return p1.Dot(p2); }

	friend TemplateVector Cross( const TemplateVector & p1, const TemplateVector & p2 )
	{ return p1 * p2; }

	TemplateVector operator*(double pR) const		// * a scalar
	{
		return TemplateVector(this->x[0]*pR, this->x[1]*pR, this->x[2]*pR);
	}

	friend TemplateVector operator*(double pR, const TemplateVector & pV)
	{
		return TemplateVector ( pV[0]*pR, pV[1]*pR, pV[2]*pR );
	}

	TemplateVector& operator*=(double pR)
	{
		for (size_t i(0); i<3; i++) this->_x[i]*=pR;
		return *this;
	}

	/// Returns the squared length
	double LenSqr(void) const
	{
		return scpr (this->getCoords (), this->getCoords (), 3);
	}

	/// Returns the length
	double Length(void) const
	{ return sqrt(LenSqr()); }

	/// Projection (component of *this parallel to pV).
	/// Note: component perpendicular to pV is:  *this - Proj(pV)
	TemplateVector Proj(const TemplateVector & pV)
	{ TemplateVector out( pV * (this->Dot(pV) / pV.LenSqr()) ); return out; }

	/// Cosine of the angle between two vectors:
	double CosAng(const TemplateVector& pV)
	{ return this->Dot(pV) / (Length() * pV.Length()); }

	/// Comparison if equal
	bool operator==( const TemplateVector & pV) const
	{
		return (std::fabs(this->_x[0] - pV[0]) < sqrt(std::numeric_limits<double>::min()) &&
				std::fabs(this->_x[1] - pV[1]) < sqrt(std::numeric_limits<double>::min()) &&
				std::fabs(this->_x[2] - pV[2]) < sqrt(std::numeric_limits<double>::min()));
	}

	/// Comparison if not equal
	bool operator!=( const TemplateVector & pV) const
	{
		return (! (pV == this));
//		this->_x[0]!=pV[0] || this->_x[1]!=pV[1] || this->_x[2]!=pV[2];
	}
};

typedef TemplateVector<double> Vector;

}

#endif // VECTOR3_H
