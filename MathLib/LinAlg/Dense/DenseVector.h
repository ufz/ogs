/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the DenseVector class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DENSEVECTOR_H_
#define DENSEVECTOR_H_

#include <vector>
#include <valarray>
#include <fstream>
#include <iterator>

namespace MathLib
{

/**
 * Dense vector class
 */
template <typename T>
class DenseVector : public std::valarray<T>
{
public:
	typedef T FP_T;

public:
	using std::valarray<T>::operator=;
	using std::valarray<T>::operator+=;
	using std::valarray<T>::operator-=;
	using std::valarray<T>::operator[];

	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 */
	explicit DenseVector(std::size_t nrows=0)
	: std::valarray<T>(nrows)
	{}

	/// return a start index of the active data range
	std::size_t getRangeBegin() const { return 0;}

	/// return an end index of the active data range
	std::size_t getRangeEnd() const { return this->size(); }

	/// get entry
	double get(std::size_t i) const { return (*this)[i]; }

	/// set a value to entry
	void set(std::size_t i, double v) { (*this)[i] = v; }

	/// add a value to entry
	void add(std::size_t i, double v) { (*this)[i] += v; }

	/**
	 * add a sub vector
	 * @param pos       positions of each sub-vector entry in this vector
	 * @param sub_vec   sub-vector
	 */
	template<class T_SUBVEC>
	void add(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec)
	{
		for (std::size_t i=0; i<pos.size(); ++i) {
			this->add(pos[i], sub_vec[i]);
		}
	}

	/**
	 * writes the matrix entries into a file
	 * @param filename output file name
	 */
	void write (const std::string &filename) const
	{
		std::ofstream os(filename);
		os << *this;
		os.close();
	}

    /// vector operation: add
    void operator+= (const DenseVector<T>& v)
    {
        *this += static_cast<std::valarray<T> >(v);
    }

    /// vector operation: subtract
    void operator-= (const DenseVector<T>& v)
    {
        *this -= static_cast<std::valarray<T> >(v);
    }

};

/**
 * writes a vector content into the output stream
 * @param os the output stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, DenseVector<T> const & v)
{
	std::copy(std::begin(v), std::end(v), std::ostream_iterator<T>(os, "\n"));
	return os;
}

}


#endif /* DENSEVECTOR_H_ */
