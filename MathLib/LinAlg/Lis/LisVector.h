/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISVECTOR_H_
#define LISVECTOR_H_

#include <iostream>

#include "lis.h"

#include "../VectorBase.h"

namespace MathLib
{

/**
 * \brief Lis vector wrapper class
 */
class LisVector : public VectorBase<double>
{
public:
    /**
     * Create a linear system using Lis
     *
     * @param length    System dimension
     */
    explicit LisVector(unsigned length);

    /**
     *
     */
    virtual ~LisVector();

    /// return a start index of the active data range
    virtual unsigned getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    virtual unsigned getRangeEnd() const { return _n_rows; }

    /// set all values in this vector
    virtual LisVector& operator= (double v)
    {
        lis_vector_set_all(v, _vec);
        return *this;
    }

    /// access entry
    virtual double operator[] (unsigned rowId) const { return get(rowId); }

    /// get entry
    virtual double get(unsigned rowId) const
    {
        double v;
        lis_vector_get_value(_vec, rowId, &v);
        return v;
    }

    /// set entry
    virtual void set(unsigned rowId, double v)
    {
        lis_vector_set_value(LIS_INS_VALUE, rowId, v, _vec);
    }

    /// add entry
    virtual void add(unsigned rowId, double v)
    {
        lis_vector_set_value(LIS_ADD_VALUE, rowId, v, _vec);
    }

    /// printout this equation for debugging
    virtual void printout(std::ostream &os=std::cout) const;

    /// return a raw Lis vector object
    LIS_VECTOR& getRawVector() {return _vec; };

private:
    LIS_VECTOR _vec;
};


} // MathLib

#endif //LISVECTOR_H_

