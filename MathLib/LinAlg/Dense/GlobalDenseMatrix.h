/**
 * @file GlobalDenseMatrix.h
 * @author Thomas Fischer
 * @date Jun 4, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef GLOBALDENSEMATRIX_H_
#define GLOBALDENSEMATRIX_H_

#include <algorithm>
#include <vector>

#include "DenseMatrix.h"

namespace MathLib
{

template<typename FP_TYPE, typename IDX_TYPE = std::size_t>
class GlobalDenseMatrix: public DenseMatrix<FP_TYPE, IDX_TYPE>
{
public:
	GlobalDenseMatrix (IDX_TYPE rows, IDX_TYPE cols) :
		DenseMatrix<FP_TYPE, IDX_TYPE>(rows, cols)
	{}

	GlobalDenseMatrix (IDX_TYPE rows, IDX_TYPE cols, const FP_TYPE& val) :
		DenseMatrix<FP_TYPE, IDX_TYPE>(rows, cols, val)
	{}
	GlobalDenseMatrix (const GlobalDenseMatrix &src) :
		DenseMatrix<FP_TYPE, IDX_TYPE>(src.getNRows(), src.getNCols())
	{
		std::copy(src._data, src._data+this->_n_rows*this->_n_cols, this->_data);
	}
	virtual ~GlobalDenseMatrix() {};

	/**
	 * Method setZero() set all matrix entries to zero.
	 */
	virtual void setZero()
	{
		std::fill(this->_data, this->_data+this->_n_rows*this->_n_cols, static_cast<FP_TYPE>(0));
	}

	/**
	 * Set the value of the matrix entry (row,col) to val.
	 * @param row The index of the row of the matrix.
	 * @param col The index of the column of the matrix.
	 * @param val The value that shoud be set.
	 * @return False if row index or column index are to large, else true.
	 */
	virtual bool setValue(IDX_TYPE row, IDX_TYPE col, FP_TYPE val)
	{
		if (row >= this->_n_rows || col >= this->_n_cols)
			return false;
		this->operator()(row,col) = val;
		return true;
	}

	/**
	 * Method adds a value to the entry at position (row,col).
	 * @param row The index of the row of the matrix.
	 * @param col The index of the column of the matrix.
	 * @param val The value that shoud be added.
	 * @return False if row index or column index are to large, else true.
	 */
	virtual bool addValue(IDX_TYPE row, IDX_TYPE col, FP_TYPE val)
	{
		if (row >= this->_n_rows || col >= this->_n_cols)
			return false;
		this->operator()(row,col) += val;
		return true;
	}

	template <class T_DENSE_MATRIX>
	void addSubMatrix(std::vector<IDX_TYPE> const& row_pos,
			std::vector<IDX_TYPE> const& col_pos, const T_DENSE_MATRIX &sub_matrix,
			FP_TYPE fkt = static_cast<FP_TYPE>(1.0))
	{
		if (row_pos.size() != sub_matrix.getNRows() || col_pos.size() != sub_matrix.getNCols())
			return;

		const std::size_t n_rows = row_pos.size();
		const std::size_t n_cols = col_pos.size();
		for (std::size_t i = 0; i < n_rows; i++) {
			const IDX_TYPE row = row_pos[i];
			for (std::size_t j = 0; j < n_cols; j++) {
				const IDX_TYPE col = col_pos[j];
				addValue(row, col, fkt * sub_matrix(i, j));
			}
		}
	}
};

} // end namespace MathLib

#endif /* GLOBALDENSEMATRIX_H_ */
