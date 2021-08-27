/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#ifndef NDEBUG
#include <fstream>
#include <string>
#endif

#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace MathLib
{
/// Global vector based on Eigen vector
class EigenVector final
{
public:
    using RawVectorType = Eigen::VectorXd;

    // The Index type of the Eigen::VectorXd class differs from the
    // Eigen::SparseMatrix<double> index type. Maybe an Eigen::SparseVector is a
    // more appropriate RawVectorType for the global vectors.
    using IndexType = Eigen::SparseMatrix<double>::Index;

    // TODO: preliminary
    EigenVector() = default;

    /// Constructor for initialization of the number of rows
    /// @param length number of rows
    explicit EigenVector(IndexType length) : vec_(length) {}

    /// copy constructor
    EigenVector(EigenVector const& src) = default;

    /// return a vector length
    IndexType size() const { return static_cast<IndexType>(vec_.size()); }

    /// return a start index of the active data range
    static constexpr IndexType getRangeBegin() { return 0; }

    /// return an end index of the active data range
    IndexType getRangeEnd() const { return size(); }

    // TODO preliminary
    void setZero() { vec_.setZero(); }

    /// access entry
    double const& operator[](IndexType rowId) const { return vec_[rowId]; }
    double& operator[](IndexType rowId) { return vec_[rowId]; }

    /// get entry
    double get(IndexType rowId) const { return vec_[rowId]; }

    /// get entries
    std::vector<double> get(std::vector<IndexType> const& indices) const
    {
        std::vector<double> local_x;
        local_x.reserve(indices.size());

        transform(cbegin(indices), cend(indices), back_inserter(local_x),
                  [&](auto const i) { return vec_[i]; });

        return local_x;
    }

    /// set entry
    void set(IndexType rowId, double v) { vec_[rowId] = v; }

    /// add entry
    void add(IndexType rowId, double v) { vec_[rowId] += v; }

    /// add entries
    template <class T_SUBVEC>
    void add(const std::vector<IndexType>& pos, const T_SUBVEC& sub_vec)
    {
        auto const length = pos.size();
        for (std::size_t i = 0; i < length; ++i)
        {
            add(pos[i], sub_vec[i]);
        }
    }

    /// Copy vector values.
    void copyValues(std::vector<double>& u) const
    {
        assert(u.size() == (std::size_t)vec_.size());
        copy_n(vec_.data(), vec_.size(), u.begin());
    }

#ifndef NDEBUG
    /// printout this equation for debugging
    void write(const std::string& filename) const
    {
        std::ofstream os(filename);
        os << vec_;
    }
#endif

    /// return a raw Eigen vector object
    RawVectorType& getRawVector() { return vec_; }

    /// return a raw Eigen vector object
    const RawVectorType& getRawVector() const { return vec_; }

private:
    RawVectorType vec_;
};

}  // namespace MathLib
