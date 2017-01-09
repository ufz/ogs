/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

#include <lis.h>

namespace MathLib
{
/**
 * \brief Lis vector wrapper class
 */
class LisVector
{
public:
    using IndexType = LIS_INT;
public:
    /**
     * Constructor for initialization of the number of rows
     * @param length number of rows
     */
    explicit LisVector(std::size_t length);

    /**
     * Constructor using the given raw data
     * @param length the length of the vector
     * @param data   the raw data
     */
    LisVector(std::size_t length, double* data);

    /// copy constructor
    LisVector(LisVector const& src);

    /**
     *
     */
    virtual ~LisVector();

    /// return a vector length
    std::size_t size() const;

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return 0; }
    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return this->size(); }

    // TODO preliminary
    void setZero() { lis_vector_set_all(0.0, _vec); }

    /// access entry
    double operator[](IndexType rowId) const { return get(rowId); }
    /// get entry
    double get(IndexType rowId) const
    {
        double v = .0;
        lis_vector_get_value(_vec, rowId, &v);
        return v;
    }

    /// set entry
    void set(IndexType rowId, double v)
    {
        lis_vector_set_value(LIS_INS_VALUE, rowId, v, _vec);
    }

    /// add entry
    void add(IndexType rowId, double v)
    {
        lis_vector_set_value(LIS_ADD_VALUE, rowId, v, _vec);
    }

    /// printout this equation for debugging
    void write(const std::string& filename) const;

    /// return a raw Lis vector object
    LIS_VECTOR& getRawVector() { return _vec; }

    ///
    template <class T_SUBVEC>
    void add(const std::vector<IndexType>& pos, const T_SUBVEC& sub_vec)
    {
        for (std::size_t i = 0; i < pos.size(); ++i)
        {
            this->add(pos[i], sub_vec[i]);
        }
    }

    /// Copy vector values.
    void copyValues(std::vector<double>& u) const
    {
        assert(u.size() == size());
        lis_vector_get_values(_vec, 0, size(), u.data());
    }

private:
    LIS_VECTOR _vec;
};

}  // MathLib
