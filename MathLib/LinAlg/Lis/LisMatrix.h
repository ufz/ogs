/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisMatrix class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISMATRIX_H_
#define LISMATRIX_H_

#include <string>
#include <vector>

#include <lis.h>

#include "MathLib/LinAlg/RowColumnIndices.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"

#include "LisOption.h"
#include "LisCheck.h"

namespace MathLib
{
// Forward declarations.
class LisVector;
class LisMatrix;

template <typename SPARSITY_PATTERN>
struct SetMatrixSparsity<LisMatrix, SPARSITY_PATTERN>;

/**
 * \brief LisMatrix is a wrapper class for matrix types of the
 * linear iterative solvers library.
 *
 * LisMatrix only supports square matrices, i.e. the number of
 * rows have to be equal to the number of columns.
 */
class LisMatrix
{
public:
    using IndexType = LIS_INT;
public:
    /**
     * constructor
     * @param n_rows the number of rows (that is equal to the number of columns)
     * @param mat_type default 1 CRS
     */
    LisMatrix(std::size_t n_rows, LisOption::MatrixType mat_type = LisOption::MatrixType::CRS);

    /**
     * constructor using raw CRS data
     *
     * Note that the given CRS data will not be deleted in this class.
     * @param n_rows   the number of rows
     * @param nonzero  the number of non zero entries
     * @param row_ptr  array of row pointer indexes
     * @param col_idx  array of column indexes
     * @param data     the non-zero entry values
     */
    LisMatrix(std::size_t n_rows, int nonzero, int* row_ptr, int* col_idx, double* data);

    /**
     *
     */
    virtual ~LisMatrix();

    /// return the number of rows
    std::size_t getNRows() const { return _n_rows; }

    /// return the number of columns
    std::size_t getNCols() const { return getNRows(); }

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return _is; }

    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return _ie; }

    /// reset this matrix with keeping its original dimension
    void setZero();

    /// set entry
    int setValue(std::size_t rowId, std::size_t colId, double v);

    /// add value
    int add(std::size_t rowId, std::size_t colId, double v);

    /// printout this equation for debugging
    void write(const std::string &filename) const;

    /// get a maximum value in diagonal entries
    double getMaxDiagCoeff();

    /// return a raw Lis matrix object
    LIS_MATRIX& getRawMatrix() { return _AA; }

    /// y = mat * x
    void multiply(const LisVector &x, LisVector &y) const;

    /// Add sub-matrix at positions \c row_pos and same column positions as the
    /// given row positions.
    template<class T_DENSE_MATRIX>
    void add(std::vector<std::size_t> const& row_pos,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(row_pos, row_pos, sub_matrix, fkt);
    }

    /// Add sub-matrix at positions given by \c indices.
    template<class T_DENSE_MATRIX>
    void add(RowColumnIndices<std::size_t> const& indices,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(indices.rows, indices.columns, sub_matrix, fkt);
    }

    ///
    template <class T_DENSE_MATRIX>
    void add(std::vector<std::size_t> const& row_pos,
            std::vector<std::size_t> const& col_pos, const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0);

    /// get this matrix type
    LisOption::MatrixType getMatrixType() const { return _mat_type; }

    /// return if this matrix is already assembled or not
    bool isAssembled() const { return _is_assembled; }

private:
    std::size_t const _n_rows;
    LisOption::MatrixType const _mat_type;
    LIS_MATRIX _AA;
    LIS_VECTOR _diag;
    bool _is_assembled;
    LIS_INT _is;	///< location where the partial matrix _AA starts in global matrix.
    LIS_INT _ie;	///< location where the partial matrix _AA ends in global matrix.
    bool _use_external_arrays;

    // friend function
    friend bool finalizeMatrixAssembly(LisMatrix &mat);

    template <typename MATRIX, typename SPARSITY_PATTERN>
    friend struct SetMatrixSparsity;
};

template<class T_DENSE_MATRIX>
void
LisMatrix::add(std::vector<std::size_t> const& row_pos, std::vector<std::size_t> const& col_pos,
        const T_DENSE_MATRIX &sub_matrix, double fkt)
{
    const std::size_t n_rows = row_pos.size();
    const std::size_t n_cols = col_pos.size();
    for (std::size_t i = 0; i < n_rows; i++) {
        const std::size_t row = row_pos[i];
        for (std::size_t j = 0; j < n_cols; j++) {
            const std::size_t col = col_pos[j];
            add(row, col, fkt * sub_matrix(i, j));
        }
    }
};

/// finish assembly to make this matrix be ready for use
bool finalizeMatrixAssembly(LisMatrix &mat);

/// Sets the sparsity pattern of the underlying LisMatrix.
template <typename SPARSITY_PATTERN>
struct SetMatrixSparsity<LisMatrix, SPARSITY_PATTERN>
{

void operator()(LisMatrix &matrix, SPARSITY_PATTERN const& sparsity_pattern)
{
    std::size_t n_rows = matrix.getNRows();
    std::vector<int> row_sizes;
    row_sizes.reserve(n_rows);

    // LIS needs 1 more entry, otherewise it starts reallocating arrays.
    for (std::size_t i = 0; i < n_rows; i++)
        row_sizes.push_back(sparsity_pattern.getNodeDegree(i) + 1);

    int ierr = lis_matrix_malloc(matrix._AA, 0, row_sizes.data());
    checkLisError(ierr);
}
};


} // MathLib

#endif //LISMATRIX_H_

