/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <fstream>
#include <string>

#include <Eigen/Sparse>

#include "MathLib/LinAlg/RowColumnIndices.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"
#include "EigenVector.h"

namespace MathLib
{

/**
 * Global matrix based on Eigen sparse matrix
 *
 * The matrix will be dynamically allocated during construction.
 */
class EigenMatrix final
{
public:
    using RawMatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    using IndexType = RawMatrixType::Index;

    // TODO The matrix constructor should take num_rows and num_cols as arguments
    //      that is left for a later refactoring.
    /**
     * constructor
     * @param n the number of rows (that is equal to the number of columns)
     * @param n_nonzero_columns the number of non-zero columns used for preallocation
     */
    explicit EigenMatrix(std::size_t n, std::size_t n_nonzero_columns = 0) :_mat(n, n)
    {
        if (n_nonzero_columns > 0)
            _mat.reserve(Eigen::VectorXi::Constant(n, n_nonzero_columns));
    }

    /// return the number of rows
    std::size_t getNumberOfRows() const { return _mat.rows(); }

    /// return the number of columns
    std::size_t getNumberOfColumns() const { return _mat.cols(); }

    /// return a start index of the active data range
    std::size_t getRangeBegin() const  { return 0; }

    /// return an end index of the active data range
    std::size_t getRangeEnd() const  { return getNumberOfRows(); }

    /// reset data entries to zero.
    void setZero()
    {
        auto const N = _mat.nonZeros();
        for (auto i = decltype(N){0}; i<N; i++)
            _mat.valuePtr()[i] = 0;
        // don't use _mat.setZero(). it makes a matrix uncompressed
    }

    /// set a value to the given entry. If the entry doesn't exist, this class
    /// dynamically allocates it.
    int setValue(IndexType row, IndexType col, double val)
    {
        assert(row < (IndexType) getNumberOfRows() && col < (IndexType) getNumberOfColumns());
        _mat.coeffRef(row, col) = val;
        return 0;
    }

    /// add a value to the given entry. If the entry doesn't exist, the value is
    /// inserted.
    int add(IndexType row, IndexType col, double val)
    {
        _mat.coeffRef(row, col) += val;
        return 0;
    }

    /// Add sub-matrix at positions \c row_pos and same column positions as the
    /// given row positions. If the entry doesn't exist, the value is inserted.
    template<class T_DENSE_MATRIX>
    void add(std::vector<IndexType> const& row_pos,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(row_pos, row_pos, sub_matrix, fkt);
    }

    /// Add sub-matrix at positions given by \c indices. If the entry doesn't exist,
    /// this class inserts the value.
    template<class T_DENSE_MATRIX>
    void add(RowColumnIndices<IndexType> const& indices,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(indices.rows, indices.columns, sub_matrix, fkt);
    }

    /// Add sub-matrix at positions \c row_pos and \c col_pos. If the entries doesn't
    /// exist in the matrix, the values are inserted.
    /// @param row_pos     a vector of row position indices. The vector size should
    ///                    equal to the number of rows in the given sub-matrix.
    /// @param col_pos     a vector of column position indices. The vector size should
    ///                    equal to the number of columns in the given sub-matrix.
    /// @param sub_matrix  a sub-matrix to be added
    /// @param fkt         a scaling factor applied to all entries in the sub-matrix
    template <class T_DENSE_MATRIX>
    void add(std::vector<IndexType> const& row_pos,
            std::vector<IndexType> const& col_pos, const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0);

    /// get value. This function returns zero if the element doesn't exist.
    double get(IndexType row, IndexType col) const
    {
        return _mat.coeff(row, col);
    }

    // TODO This method is currently used nowhere.
    //      I disabled it in order to simplify the interface of this class,
    //      especially since our other matrix classes lack this operator.
    /*
    /// get value. This function returns zero if the element doesn't exist.
    double operator() (IndexType row, IndexType col) const
    {
        return get(row, col);
    }
    */

    /// get a maximum value in diagonal entries
    double getMaxDiagCoeff() const
    {
        return _mat.diagonal().maxCoeff();
    }

    /// return always true, i.e. the matrix is always ready for use
    bool isAssembled() const { return true; }

    /// printout this matrix for debugging
    void write(const std::string &filename) const
    {
        std::ofstream of(filename);
        if (of)
            write(of);
    }

    /// printout this matrix for debugging
    void write(std::ostream &os) const
    {
        for (int k=0; k<_mat.outerSize(); ++k)
          for (RawMatrixType::InnerIterator it(_mat,k); it; ++it)
              os << it.row() << " " << it.col() << ": " << it.value() << "\n";
        os << std::endl;
    }

    RawMatrixType& getRawMatrix() { return _mat; }
    const RawMatrixType& getRawMatrix() const { return _mat; }

protected:
    RawMatrixType _mat;
};

template <class T_DENSE_MATRIX>
void EigenMatrix::add(std::vector<IndexType> const& row_pos,
                      std::vector<IndexType> const& col_pos,
                      const T_DENSE_MATRIX& sub_matrix, double fkt)
{
    auto const n_rows = row_pos.size();
    auto const n_cols = col_pos.size();
    for (auto i = decltype(n_rows){0}; i < n_rows; i++) {
        auto const row = row_pos[i];
        for (auto j = decltype(n_cols){0}; j < n_cols; j++) {
            auto const col = col_pos[j];
            add(row, col, fkt * sub_matrix(i, j));
        }
    }
};

/// Sets the sparsity pattern of the underlying EigenMatrix.
template <typename SPARSITY_PATTERN>
struct SetMatrixSparsity<EigenMatrix, SPARSITY_PATTERN>
{

/// \note This operator relies on row-major storage order of the underlying
/// eigen matrix i.e. of the RawMatrixType.
void operator()(EigenMatrix &matrix, SPARSITY_PATTERN const& sparsity_pattern)
{
    static_assert(EigenMatrix::RawMatrixType::IsRowMajor,
                  "Set matrix sparsity relies on the EigenMatrix to be in "
                  "row-major storage order.");

    assert(matrix.getNumberOfRows() == sparsity_pattern.size());

    matrix.getRawMatrix().reserve(sparsity_pattern);
}
};

} // end namespace MathLib
