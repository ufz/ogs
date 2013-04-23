/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the Vector class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <valarray>

#include "../VectorBase.h"


namespace MathLib
{

/**
 * class VectorBase is the basis for all vector classes
 */
template <typename T>
class Vector : public VectorBase<T>
{
    typedef Vector<T> _vector_type;
public:
	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 * @return
	 */
    Vector(unsigned nrows=0) :
        VectorBase<T>(nrows), _data(nrows)
	{}

	/**
	 * copy constructor.
	 * @param original the object that is copied
	 * @return
	 */
    Vector (Vector const& original) :
        VectorBase<T>(original), _data(original._data)
	{}

	/**
	 * destructor of the class.
	 * @return
	 */
	virtual ~Vector() {};

    /// get entry
    virtual double get(unsigned i) const { return _data[i]; };

    /// set a value to entry
    virtual void set(unsigned i, double v) { _data[i] = v; };

    /// add a value to entry
    virtual void add(unsigned i, double v) { _data[i] += v; };

    /// access entry
    T& operator[] (unsigned i) { return _data[i]; }

    /// access entry
    virtual T operator[] (unsigned i) const { return _data[i]; };

    /// vector operation: set data
    _vector_type& operator= (const _vector_type &src)
    {
        _data = src._data;
        return *this;
    }

    /// vector operation: add
    void operator+= (const _vector_type& v) { _data += v._data; };

    /// vector operation: subtract
    void operator-= (const _vector_type& v) { _data -= v._data; };

    /// set all values in this vector
    _vector_type& operator= (T val)
    {
        _data = val;
        return *this;
    }

    /// return pointer to raw data
    T* getData() { return &_data[0]; }

    /// return pointer to raw data
    const T* getData() const { return &_data[0]; }

    /**
     * writes the matrix entries into the output stream
     * @param out the output stream
     */
    void write (std::ostream& out) const
    {
        for (std::size_t i = 0; i < _data.size(); i++) {
            out << _data[i] << "\n";
        }
    }

protected:
    std::valarray<T> _data;
};

}


#endif /* VECTOR_H_ */
