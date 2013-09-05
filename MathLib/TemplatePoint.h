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
#include <array>
#include <algorithm>
#include <iterator>
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
template <typename T, std::size_t DIM = 3> class TemplatePoint
{
public:
	/** default constructor */
	TemplatePoint();

	/** constructor - constructs a TemplatePoint object
	 *
	 * @param x std::array containing the coordinates of the point
	 */
	TemplatePoint(std::array<T,DIM> const& x);

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
	const T* getCoords () const { return &_x[0]; }

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
		for (std::size_t k(0); k<DIM; k++)
			is >> _x[k];
	}

protected:
	std::array<T, DIM> _x;
};

template <typename T, std::size_t DIM>
TemplatePoint<T,DIM>::TemplatePoint()
{}

template <typename T, std::size_t DIM>
TemplatePoint<T,DIM>::TemplatePoint(std::array<T,DIM> const& x) :
	_x(x)
{}

template <typename T, std::size_t DIM>
TemplatePoint<T, DIM>::TemplatePoint (T const* x)
{
	for (std::size_t k(0); k < DIM; k++)
		_x[k] = x[k];
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
