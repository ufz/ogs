/*
 * CRSMatrix.h
 *
 *  Created on: Sep 20, 2011
 *      Author: TF
 */

#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

// Base
#include "swap.h"

// MathLib
#include "SparseMatrixBase.h"
#include "sparse.h"
#include "amuxCRS.h"
#include "../Preconditioner/generateDiagPrecond.h"

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE>
class CRSMatrix: public SparseMatrixBase<FP_TYPE, IDX_TYPE>
{
public:
	CRSMatrix(std::string const &fname) :
		SparseMatrixBase<FP_TYPE, IDX_TYPE>(),
		_row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{
		std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
		if (in) {
			CS_read(in, SparseMatrixBase<FP_TYPE, IDX_TYPE>::_n_rows, _row_ptr, _col_idx, _data);
			SparseMatrixBase<FP_TYPE, IDX_TYPE>::_n_cols = SparseMatrixBase<FP_TYPE, IDX_TYPE>::_n_rows;
			in.close();
		} else {
			std::cout << "cannot open " << fname << std::endl;
		}
	}

	CRSMatrix(IDX_TYPE n, IDX_TYPE *iA, IDX_TYPE *jA, FP_TYPE* A) :
		SparseMatrixBase<FP_TYPE, IDX_TYPE>(n,n),
		_row_ptr(iA), _col_idx(jA), _data(A)
	{}

	CRSMatrix(IDX_TYPE n1) :
		SparseMatrixBase<FP_TYPE, IDX_TYPE>(n1, n1),
		_row_ptr(NULL), _col_idx(NULL), _data(NULL)
	{}

	virtual ~CRSMatrix()
	{
		delete [] _row_ptr;
		delete [] _col_idx;
		delete [] _data;
	}

	virtual void amux(FP_TYPE d, FP_TYPE const * const __restrict__ x, FP_TYPE * __restrict__ y) const
	{
		amuxCRS<FP_TYPE, IDX_TYPE>(d, this->getNRows(), _row_ptr, _col_idx, _data, x, y);
	}

    virtual void precondApply(FP_TYPE* /*x*/) const
    {}

    /**
     * get the number of non-zero entries
     * @return number of non-zero entries
     */
    IDX_TYPE getNNZ() const { return _row_ptr[MatrixBase::_n_rows]; }

    /**
     * This method inserts/overwrites a non-zero matrix entry.
	 * Precondition: the entry have to be in the sparsity pattern!
     * @param row the row number
     * @param col the column number
     * @param val the value that should be set at pos row,col
     * @return a value > 0, if the entry is not contained in the sparsity pattern
     */
	int setValue(IDX_TYPE row, IDX_TYPE col, FP_TYPE val)
	{
		assert(0 <= row && row < MatrixBase::_n_rows);

		// linear search - for matrices with many entries per row binary search is much faster
		const IDX_TYPE idx_end (_row_ptr[row+1]);
		IDX_TYPE j(_row_ptr[row]), k;

		while (j<idx_end && (k=_col_idx[j]) <= col) {
			if (k == col) {
				_data[j] = val;
				return 0;
			}
			j++;
		}
		return 1;
	}

    /**
     * This method adds value val to an existing matrix entry at position row,col.
	 * Precondition: the entry have to be in the sparsity pattern!
     * @param row the row number
     * @param col the column number
     * @param val the value that should be set at pos row,col
     * @return a value > 0, if the entry is not contained in the sparsity pattern
     */
	int addValue(IDX_TYPE row, IDX_TYPE col, FP_TYPE val)
	{
		assert(0 <= row && row < MatrixBase::_n_rows);

		// linear search - for matrices with many entries per row binary search is much faster
		const IDX_TYPE idx_end (_row_ptr[row+1]);
		IDX_TYPE j(_row_ptr[row]), k;

		while (j<idx_end && (k=_col_idx[j]) <= col) {
			if (k == col) {
				#pragma omp atomic
				_data[j] += val;
				return 0;
			}
			j++;
		}
		return 1;
	}

    /**
     * This is an access operator to a non-zero matrix entry. If the value of
     * a non-existing matrix entry is requested it will be 0.0 returned.
     * @param row the row number
     * @param col the column number
     * @return The corresponding matrix entry or 0.0.
     */
	double getValue(IDX_TYPE row, IDX_TYPE col)
	{
		assert(0 <= row && row < MatrixBase::_n_rows);

		// linear search - for matrices with many entries per row binary search is much faster
		const IDX_TYPE idx_end (_row_ptr[row+1]);
		IDX_TYPE j(_row_ptr[row]), k;

		while (j<idx_end && (k=_col_idx[j]) <= col) {
			if (k == col) {
				return _data[j];
			}
			j++;
		}
		return 0.0;
	}

    /**
     * This is the constant access operator to a non-zero matrix entry.
     * Precondition: the entries have to be in the sparsity pattern!
     * @param row the row number
     * @param col the column number
     * @return The corresponding matrix entry.
     */
    FP_TYPE operator() (IDX_TYPE row, IDX_TYPE col) const
    {
    	assert(0 <= row && row < MatrixBase::_n_rows);

    	// linear search - for matrices with many entries per row binary search is much faster
    	const IDX_TYPE idx_end (_row_ptr[row+1]);
    	IDX_TYPE j(_row_ptr[row]), k;

    	while (j<idx_end && (k=_col_idx[j]) <= col) {
    		if (k == col) {
    			return _data[j];
    		}
    		j++;
    	}
    	return 0.0;
    }

	/**
	 * get const access to the row pointer array of CRS matrix
	 * @return the index array _row_ptr
	 */
	IDX_TYPE const* getRowPtrArray() const { return _row_ptr; }

	/**
	 * get const access to the column index array of CRS matrix
	 * @return the index array _col_idx
	 */
	IDX_TYPE const* getColIdxArray ()const { return _col_idx; }

	/**
	 * get the matrix entries within an array of CRS matrix
	 * @return
	 */
	FP_TYPE const* getEntryArray() const { return _data; }

	/**
	 * erase rows and columns from sparse matrix
	 * @param n_rows_cols number of rows / columns to remove
	 * @param rows_cols sorted list of rows/columns that should be removed
	 */
	void eraseEntries(IDX_TYPE n_rows_cols, IDX_TYPE const* const rows_cols)
	{
		//*** remove the rows
		removeRows(n_rows_cols, rows_cols);
		//*** transpose
		transpose();
		//*** remove columns in original means removing rows in the transposed
		removeRows(n_rows_cols, rows_cols);
		//*** transpose again
		transpose();
	}

	/**
	 * get the j-th column of the sparse matrix
	 * @param j the column number that should be returned
	 * @param column_entries the column entries (have to be allocated
	 */
	void getColumn(IDX_TYPE j, FP_TYPE* column_entries) const
	{
		for (IDX_TYPE k(0); k<MatrixBase::_n_rows; k++) {
			const IDX_TYPE end_row(_row_ptr[k+1]);
			IDX_TYPE i(_row_ptr[k+1]);
			while (i<end_row && _col_idx[i] != j) {
				i++;
			}
			if (i==end_row) {
				column_entries[k] = 0.0;
			} else {
				column_entries[k] = _data[i];
			}
		}
	}

protected:
	void removeRows (IDX_TYPE n_rows_cols, IDX_TYPE const*const rows)
	{
		//*** determine the number of new rows and the number of entries without the rows
		const IDX_TYPE n_new_rows(MatrixBase::_n_rows - n_rows_cols);
		IDX_TYPE *row_ptr_new(new IDX_TYPE[n_new_rows+1]);
		row_ptr_new[0] = 0;
		IDX_TYPE row_cnt (1), erase_row_cnt(0);
		for (unsigned k(0); k<MatrixBase::_n_rows; k++) {
			if (erase_row_cnt < n_rows_cols) {
				if (k != rows[erase_row_cnt]) {
					row_ptr_new[row_cnt] = _row_ptr[k+1] - _row_ptr[k];
					row_cnt++;
				} else {
					erase_row_cnt++;
				}
			} else {
				row_ptr_new[row_cnt] = _row_ptr[k+1] - _row_ptr[k];
				row_cnt++;
			}
		}

		//*** sum up the entries
		for (IDX_TYPE k(0); k<n_new_rows; k++) {
			row_ptr_new[k+1] = row_ptr_new[k+1] + row_ptr_new[k];
		}

		//*** create new memory for col_idx and data
		IDX_TYPE nnz_new(row_ptr_new[n_new_rows]);
		IDX_TYPE *col_idx_new (new IDX_TYPE[nnz_new]);
		FP_TYPE *data_new (new FP_TYPE[nnz_new]);

		//*** copy the entries
		// initialization
		IDX_TYPE *row_ptr_new_tmp(new IDX_TYPE[n_new_rows+1]);
		for (unsigned k(0); k<=n_new_rows; k++) {
			row_ptr_new_tmp[k] = row_ptr_new[k];
		}
		erase_row_cnt = 0;
		row_cnt = 0;
		// copy column index and data entries
		for (IDX_TYPE k(0); k<MatrixBase::_n_rows; k++) {
			if (erase_row_cnt < n_rows_cols) {
				if (k != rows[erase_row_cnt]) {
					const IDX_TYPE end (_row_ptr[k+1]);
					// walk through row
					for (IDX_TYPE j(_row_ptr[k]); j<end; j++) {
						col_idx_new[row_ptr_new_tmp[row_cnt]] = _col_idx[j];
						data_new[row_ptr_new_tmp[row_cnt]] = _data[j];
						row_ptr_new_tmp[row_cnt]++;
					}
					row_cnt++;
				} else {
					erase_row_cnt++;
				}
			} else {
				const IDX_TYPE end (_row_ptr[k+1]);
				// walk through row
				for (IDX_TYPE j(_row_ptr[k]); j<end; j++) {
					col_idx_new[row_ptr_new_tmp[row_cnt]] = _col_idx[j];
					data_new[row_ptr_new_tmp[row_cnt]] = _data[j];
					row_ptr_new_tmp[row_cnt]++;
				}
				row_cnt++;
			}
		}

		MatrixBase::_n_rows -= n_rows_cols;
		BaseLib::swap (row_ptr_new, _row_ptr);
		BaseLib::swap (col_idx_new, _col_idx);
		BaseLib::swap (data_new, _data);

		delete [] row_ptr_new_tmp;
		delete [] row_ptr_new;
		delete [] col_idx_new;
		delete [] data_new;
	}

	void transpose ()
	{
		// create a helper array row_ptr_nnz
		IDX_TYPE *row_ptr_nnz(new IDX_TYPE[MatrixBase::_n_cols+1]);
		for (IDX_TYPE k(0); k <= MatrixBase::_n_cols; k++) {
			row_ptr_nnz[k] = 0;
		}

		// count entries per row in the transposed matrix
		IDX_TYPE nnz(_row_ptr[MatrixBase::_n_rows]);
		for (IDX_TYPE k(0); k < nnz; k++) {
			row_ptr_nnz[_col_idx[k]]++;
		}

		// create row_ptr_trans
		IDX_TYPE *row_ptr_trans(new IDX_TYPE[MatrixBase::_n_cols + 1]);
		row_ptr_trans[0] = 0;
		for (IDX_TYPE k(0); k < MatrixBase::_n_cols; k++) {
			row_ptr_trans[k+1] = row_ptr_trans[k] + row_ptr_nnz[k];
		}

		// make a copy of row_ptr_trans
		for (IDX_TYPE k(0); k <= MatrixBase::_n_cols; k++) {
			row_ptr_nnz[k] = row_ptr_trans[k];
		}

		// create arrays col_idx_trans and data_trans
		assert(nnz == row_ptr_trans[MatrixBase::_n_cols]);
		IDX_TYPE *col_idx_trans(new IDX_TYPE[nnz]);
		FP_TYPE *data_trans(new FP_TYPE[nnz]);

		// fill arrays col_idx_trans and data_trans
		for (IDX_TYPE i(0); i < MatrixBase::_n_rows; i++) {
			const IDX_TYPE row_end(_row_ptr[i + 1]);
			for (IDX_TYPE j(_row_ptr[i]); j < row_end; j++) {
				const IDX_TYPE k(_col_idx[j]);
				col_idx_trans[row_ptr_nnz[k]] = i;
				data_trans[row_ptr_nnz[k]] = _data[j];
				row_ptr_nnz[k]++;
			}
		}

		BaseLib::swap(MatrixBase::_n_rows, MatrixBase::_n_cols);
		BaseLib::swap(row_ptr_trans, _row_ptr);
		BaseLib::swap(col_idx_trans, _col_idx);
		BaseLib::swap(data_trans, _data);

		delete[] row_ptr_nnz;
		delete[] row_ptr_trans;
		delete[] col_idx_trans;
		delete[] data_trans;
	}

#ifndef NDEBUG
	void printMat() const
	{
		for (IDX_TYPE k(0); k<MatrixBase::_n_rows; k++) {
			std::cout << k << ": " << std::flush;
			const IDX_TYPE row_end(_row_ptr[k+1]);
			for (IDX_TYPE j(_row_ptr[k]); j<row_end; j++) {
				std::cout << _col_idx[j] << " " << std::flush;
			}
			std::cout << std::endl;
		}
	}
#endif

	IDX_TYPE *_row_ptr;
	IDX_TYPE *_col_idx;
	FP_TYPE* _data;
};

} // end namespace MathLib

#endif

