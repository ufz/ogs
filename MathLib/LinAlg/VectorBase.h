/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the VectorBase class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTORBASE_H_
#define VECTORBASE_H_

#include <vector>

namespace MathLib
{

/**
 * class VectorBase is the basis for all vector classes
 */
template <typename T>
class VectorBase
{
    typedef VectorBase<T> _vector_type;
public:
	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 * @return
	 */
    explicit VectorBase(unsigned nrows=0) :
		_n_rows(nrows)
	{}

	/**
	 * copy constructor.
	 * @param original the object that is copied
	 * @return
	 */
    VectorBase (VectorBase const& original) :
		_n_rows (original._n_rows)
	{}

	/**
	 * destructor of the class.
	 * @return
	 */
	virtual ~VectorBase() {};

	/**
	 * get the number of rows
	 * @return the number of rows
	 */
	unsigned getNRows () const { return _n_rows; }

    /// return a start index of the active data range
    virtual unsigned getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    virtual unsigned getRangeEnd() const { return _n_rows; }

    /// get entry
    virtual double get(unsigned rowId) const = 0;

    /// set a value to entry
    virtual void set(unsigned rowId, double v) = 0;

    /// add a value to entry
    virtual void add(unsigned rowId, double v) = 0;

    /// access entry
    virtual T operator[] (unsigned idx) const = 0;

    /// vector operation: set data
    virtual _vector_type& operator= (const _vector_type &src);

    /// vector operation: add
    virtual void operator+= (const _vector_type& v);

    /// vector operation: subtract
    virtual void operator-= (const _vector_type& v);

    /// set all values in this vector
    virtual _vector_type& operator= (T v);

    template<class T_SUBVEC>
    void addSubVector(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec)
    {
        for (std::size_t i=0; i<pos.size(); ++i) {
            this->add(pos[i], sub_vec[i]);
        }
    }

protected:
	/**
	 * the number of rows
	 */
	unsigned _n_rows;
};

}

#include "VectorBase.tpp"

#endif /* VECTORBASE_H_ */
