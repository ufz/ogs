/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisMatrix class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    /// Matrix type
    enum class MatrixType : int
    {
        CRS = 1,
        CCS = 2,
        MSR = 3,
        DIA = 4,
        ELL = 5,
        JDS = 6,
        BSR = 7,
        BSC = 8,
        VBR = 9,
        COO = 10,
        DNS = 11
    };

    using IndexType = LIS_INT;
public:
    /**
     * constructor
     * @param n_rows the number of rows (that is equal to the number of columns)
     * @param mat_type default 1 CRS
     */
    LisMatrix(std::size_t n_rows, MatrixType mat_type = MatrixType::CRS);

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
    LisMatrix(std::size_t n_rows, int nonzero, IndexType* row_ptr, IndexType* col_idx,
              double* data);

    /**
     *
     */
    virtual ~LisMatrix();

    /// return the number of rows
    std::size_t getNumberOfRows() const { return _n_rows; }

    /// return the number of columns
    std::size_t getNumberOfColumns() const { return getNumberOfRows(); }

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return _is; }

    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return _ie; }

    /// reset this matrix with keeping its original dimension
    void setZero();

    /// set entry
    int setValue(IndexType rowId, IndexType colId, double v);

    /// add value
    int add(IndexType rowId, IndexType colId, double v);

    /// printout this equation for debugging
    void write(const std::string &filename) const;

    /// get a maximum value in diagonal entries
    double getMaxDiagCoeff();

    /// return a raw Lis matrix object
    LIS_MATRIX& getRawMatrix() { return _AA; }

    /// Add sub-matrix at positions \c row_pos and same column positions as the
    /// given row positions.
    template<class T_DENSE_MATRIX>
    void add(std::vector<IndexType> const& row_pos,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(row_pos, row_pos, sub_matrix, fkt);
    }

    /// Add sub-matrix at positions given by \c indices.
    template<class T_DENSE_MATRIX>
    void add(RowColumnIndices<IndexType> const& indices,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(indices.rows, indices.columns, sub_matrix, fkt);
    }

    template <class T_DENSE_MATRIX>
    void add(std::vector<IndexType> const& row_pos,
             std::vector<IndexType> const& col_pos,
             const T_DENSE_MATRIX& sub_matrix, double fkt = 1.0);

    /// get this matrix type
    MatrixType getMatrixType() const { return _mat_type; }

    /// return if this matrix is already assembled or not
    bool isAssembled() const { return _is_assembled; }

private:
    std::size_t const _n_rows;
    MatrixType const _mat_type;
    LIS_MATRIX _AA;
    LIS_VECTOR _diag;
    bool _is_assembled;
    IndexType _is;    ///< location where the partial matrix _AA starts in global matrix.
    IndexType _ie;    ///< location where the partial matrix _AA ends in global matrix.
    bool _use_external_arrays;

    // friend function
    friend bool finalizeMatrixAssembly(LisMatrix &mat);

    template <typename MATRIX, typename SPARSITY_PATTERN>
    friend struct SetMatrixSparsity;
};

template <class T_DENSE_MATRIX>
void LisMatrix::add(std::vector<IndexType> const& row_pos,
                    std::vector<IndexType> const& col_pos,
                    const T_DENSE_MATRIX& sub_matrix, double fkt)
{
    auto const n_rows = row_pos.size();
    auto const n_cols = col_pos.size();
    for (auto i = decltype(n_rows){0}; i < n_rows; i++)
    {
        auto const row = row_pos[i];
        for (auto j = decltype(n_cols){0}; j < n_cols; j++)
        {
            auto const col = col_pos[j];
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
    auto const n_rows = matrix.getNumberOfRows();
    std::vector<LisMatrix::IndexType> row_sizes;
    row_sizes.reserve(n_rows);

    // LIS needs 1 more entry, otherewise it starts reallocating arrays.
    for (auto i : sparsity_pattern) row_sizes.push_back(i+1);

    int ierr = lis_matrix_malloc(matrix._AA, 0, row_sizes.data());
    checkLisError(ierr);
}
};


} // MathLib

#endif //LISMATRIX_H_
