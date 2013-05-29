/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the DenseVector class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace MathLib
{

/**
 * Dense vector class
 */
template <typename T>
class DenseVector : public std::valarray<T>
{
public:
	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 * @return
	 */
	explicit DenseVector(std::size_t nrows=0)
	: std::valarray<T>(nrows)
	{}

	/**
	 * copy constructor.
	 * @param original the object that is copied
	 * @return
	 */
	DenseVector (DenseVector<T> const& original)
	: std::valarray<T>(static_cast<std::valarray<T> >(original))
	{}

	/**
	 * destructor of the class.
	 * @return
	 */
	virtual ~DenseVector() {};

	/**
	 * set a given value to all entries in this vector
	 *
	 * @param val   value
	 * @return a reference to this vector
	 */
	DenseVector<T> & operator=(const T& val)
	{
		this->std::valarray<T>::operator =(val);
		return *this;
	};

	/// return a start index of the active data range
	std::size_t getRangeBegin() const { return 0;}

	/// return an end index of the active data range
	std::size_t getRangeEnd() const { return this->size(); }

	/// get entry
	double get(std::size_t i) const { return (*this)[i]; };

	/// set a value to entry
	void set(std::size_t i, double v) { (*this)[i] = v; };

	/// add a value to entry
	void add(std::size_t i, double v) { (*this)[i] += v; };

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

	/**
	 * writes the matrix entries into a file
	 * @param output file name
	 */
	void write (const std::string &filename) const
	{
		std::ofstream os(filename);
		os << *this;
		os.close();
	}

};

/**
 * writes a vector content into the output stream
 * @param out the output stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, DenseVector<T> const & v)
{
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, "\n"));
	return os;
}

}


#endif /* DENSEVECTOR_H_ */
