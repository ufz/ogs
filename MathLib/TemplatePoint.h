/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-28
 * \brief  Definition of the TemplatePoint class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEPOINT_H_
#define TEMPLATEPOINT_H_

// STL
#include <cassert>
#include <iostream>

namespace MathLib
{
/**
 * \ingroup GeoLib
 *
 * \brief class-template for points can be instantiated by a numeric type.
 * \tparam T the coordinate type
 */
template <class T> class TemplatePoint
{
public:
	/** default constructor */
	TemplatePoint();

	/** constructor - constructs a TemplatePoint object
	   \param x1 value for the first coordinate
	   \param x2 value for the second coordinate
	   \param x3 value for the third coordinate
	 */
	TemplatePoint(T x1, T x2, T x3);

	/** constructor - constructs a TemplatePoint object
	   \param x values for three coordinates
	 */
	TemplatePoint (T const* x);

	/** virtual destructor */
	virtual ~TemplatePoint() {}

	/** \brief const access operator
	 *  The access to the point coordinates is like the access to a field. Code example:
	 * \code
	 * Point<double> point (1.0, 2.0, 3.0);
	 * double sqrNrm2 = point[0] * point[0] + point[1] * point[1] + point[2] + point[2];
	 * \endcode
	 */
	const T& operator[] (std::size_t idx) const
	{
		assert (idx <= 2);
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
	const T* getCoords () const { return _x; }

	/** write point coordinates into stream (used from operator<<)
	 * \param os a standard output stream
	 */
	virtual void write (std::ostream &os) const
	{
		os << _x[0] << " " << _x[1] << " " << _x[2] << std::flush;
	}

	/** read point coordinates into stream (used from operator>>) */
	virtual void read (std::istream &is)
	{
		is >> _x[0] >> _x[1] >> _x[2];
	}

protected:
	T _x[3];
};

template <class T> TemplatePoint<T>::TemplatePoint()
{
	_x[0] = static_cast<T>(0);
	_x[1] = static_cast<T>(0);
	_x[2] = static_cast<T>(0);
}

template <class T> TemplatePoint<T>::TemplatePoint(T x1, T x2, T x3)
{
	_x[0] = x1;
	_x[1] = x2;
	_x[2] = x3;
}

template <class T> TemplatePoint<T>::TemplatePoint (T const* x)
{
	for (std::size_t k(0); k < 3; k++)
		_x[k] = x[k];
}

/** overload the output operator for class Point */
template <class T>
std::ostream& operator<< (std::ostream &os, const TemplatePoint<T> &p)
{
	p.write (os);
	return os;
}

/** overload the input operator for class Point */
template <class T>
std::istream& operator>> (std::istream &is, TemplatePoint<T> &p)
{
	p.read (is);
	return is;
}
} // end namespace MathLib

#endif /* TEMPLATEPOINT_H_ */
