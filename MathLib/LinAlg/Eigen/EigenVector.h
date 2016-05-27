/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENVECTOR_H_
#define EIGENVECTOR_H_

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
    EigenVector() {}

    /// Constructor for initialization of the number of rows
    /// @param length number of rows
    explicit EigenVector(std::size_t length) : _vec(length) {}

    /// copy constructor
    EigenVector(EigenVector const &src) : _vec(src._vec) {}

    /// return a vector length
    std::size_t size() const { return _vec.size(); }

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return size(); }

    /// set all values in this vector
    EigenVector& operator= (double v) { _vec.setConstant(v); return *this; }

    // TODO preliminary
    void setZero() { _vec.setZero(); }

    /// set all values in this vector
    EigenVector& operator*= (double v) { _vec *= v; return *this; }

    /// access entry
    double const & operator[] (IndexType rowId) const { return _vec[rowId]; }
    double& operator[] (IndexType rowId) { return _vec[rowId]; }

    /// get entry
    double get(IndexType rowId) const
    {
        return _vec[rowId];
    }

    /// set entry
    void set(IndexType rowId, double v)
    {
        _vec[rowId] = v;
    }

    /// add entry
    void add(IndexType rowId, double v)
    {
        _vec[rowId] += v;
    }

    /// add entries
    template<class T_SUBVEC>
    void add(const std::vector<IndexType> &pos, const T_SUBVEC &sub_vec)
    {
        for (std::size_t i=0; i<pos.size(); ++i) {
            add(pos[i], sub_vec[i]);
        }
    }

    /// Copy vector values.
    void copyValues(std::vector<double>& u) const
    {
        assert(u.size() == (std::size_t) _vec.size());
        copy_n(_vec.data(), _vec.size(), u.begin());
    }

#ifndef NDEBUG
    /// printout this equation for debugging
    void write (const std::string &filename) const {
		std::ofstream os(filename);
		os.precision(16);
		os << _vec;
	}
#endif

    /// return a raw Eigen vector object
    RawVectorType& getRawVector() {return _vec; }

    /// return a raw Eigen vector object
    const RawVectorType& getRawVector() const {return _vec; }

    /// vector operation: set data
    EigenVector& operator= (const EigenVector &src) { _vec = src._vec; return *this; }

    /// vector operation: add
    void operator+= (const EigenVector& v) { _vec += v._vec; }

    /// vector operation: subtract
    void operator-= (const EigenVector& v) { _vec -= v._vec; }

private:
    RawVectorType _vec;
};

} // MathLib

#endif //EIGENVECTOR_H_

