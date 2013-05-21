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
class Vector : public std::valarray<T>
{
public:
	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 * @return
	 */
    explicit Vector(unsigned nrows=0)
    : std::valarray<T>(nrows)
	{}

	/**
	 * copy constructor.
	 * @param original the object that is copied
	 * @return
	 */
    Vector (Vector<T> const& original)
    : std::valarray<T>(static_cast<std::valarray<T> >(original))
	{}

	/**
	 * destructor of the class.
	 * @return
	 */
	virtual ~Vector() {};

	/**
	 * set a given value to all entries in this vector
	 *
	 * @param val   value
	 * @return a reference to this vector
	 */
	Vector<T> & operator=(const T& val)
	{
	    this->std::valarray<T>::operator =(val);
	    return *this;
	};

    /// return a start index of the active data range
    unsigned getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    unsigned getRangeEnd() const { return this->size(); }

    /// get entry
    double get(unsigned i) const { return (*this)[i]; };

    /// set a value to entry
    void set(unsigned i, double v) { (*this)[i] = v; };

    /// add a value to entry
    void add(unsigned i, double v) { (*this)[i] += v; };

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
     * writes the matrix entries into the output stream
     * @param out the output stream
     */
    void write (std::ostream& out) const
    {
        for (std::size_t i = 0; i < this->size(); i++) {
            out << (*this)[i] << "\n";
        }
    }

};

}


#endif /* VECTOR_H_ */
