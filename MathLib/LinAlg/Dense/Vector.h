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

#include <vector>
#include <valarray>


namespace MathLib
{

/**
 * Dense vector class
 */
template <typename T>
class Vector
{
    typedef Vector<T> _vector_type;
public:
	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 * @return
	 */
    Vector(unsigned nrows=0)
    : _data(nrows)
	{}

	/**
	 * copy constructor.
	 * @param original the object that is copied
	 * @return
	 */
    Vector (Vector const& original)
    : _data(original._data)
	{}

	/**
	 * destructor of the class.
	 * @return
	 */
	virtual ~Vector() {};

    /**
     * get the number of rows
     * @return the number of rows
     */
    unsigned getNRows() const { return _data.size(); }

    /// return a start index of the active data range
    unsigned getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    unsigned getRangeEnd() const { return getNRows(); }

    /// get entry
    double get(unsigned i) const { return _data[i]; };

    /// set a value to entry
    void set(unsigned i, double v) { _data[i] = v; };

    /// add a value to entry
    void add(unsigned i, double v) { _data[i] += v; };

    /**
     * add a sub vector
     * @param pos       positions of each sub-vector entry in this vector
     * @param sub_vec   sub-vector
     */
    template<class T_SUBVEC>
    void addSubVector(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec)
    {
        for (std::size_t i=0; i<pos.size(); ++i) {
            this->add(pos[i], sub_vec[i]);
        }
    }

    /// access entry
    T& operator[] (unsigned i) { return _data[i]; }

    /// access entry
    T operator[] (unsigned i) const { return _data[i]; };

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
