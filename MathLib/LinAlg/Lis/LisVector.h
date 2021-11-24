/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <lis.h>

#include <cassert>
#include <string>
#include <vector>

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
     */
    virtual ~LisVector();

    /// return a vector length
    std::size_t size() const;

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return 0; }
    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return this->size(); }

    // TODO preliminary
    void setZero() { lis_vector_set_all(0.0, vec_); }

    /// access entry
    double operator[](IndexType rowId) const { return get(rowId); }
    /// get entry
    double get(IndexType rowId) const
    {
        double v = .0;
        lis_vector_get_value(vec_, rowId, &v);
        return v;
    }

    /// set entry
    void set(IndexType rowId, double v)
    {
        lis_vector_set_value(LIS_INS_VALUE, rowId, v, vec_);
    }

    /// add entry
    void add(IndexType rowId, double v)
    {
        lis_vector_set_value(LIS_ADD_VALUE, rowId, v, vec_);
    }

    /// printout this equation for debugging
    void write(const std::string& filename) const;

    /// return a raw Lis vector object
    LIS_VECTOR& getRawVector() { return vec_; }

    ///
    template <class T_SUBVEC>
    void add(const std::vector<IndexType>& pos, const T_SUBVEC& sub_vec)
    {
        for (std::size_t i = 0; i < pos.size(); ++i)
        {
            this->add(pos[i], sub_vec[i]);
        }
    }

    /// Copy local entries to a vector.
    /// \param u a vector for the values of local entries. It will be resized to
    /// hold the current vector data.
    void copyValues(std::vector<double>& u) const
    {
        u.resize(size());
        lis_vector_get_values(vec_, 0, size(), u.data());
    }

private:
    LIS_VECTOR vec_;
};

}  // MathLib
