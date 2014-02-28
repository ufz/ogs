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
#include "DenseVector.h"

#include "MathLib/LinAlg/RowColumnIndices.h"

namespace MathLib
{

template<typename FP_TYPE, typename IDX_TYPE = std::size_t>
class GlobalDenseMatrix: public DenseMatrix<FP_TYPE, IDX_TYPE>
{
public:
	typedef FP_TYPE FP_T;

public:
	/// Dense square matrix constructor.
	GlobalDenseMatrix (IDX_TYPE rows);

	/// Dense rectangular matrix constructor.
	GlobalDenseMatrix (IDX_TYPE rows, IDX_TYPE cols);
	GlobalDenseMatrix (IDX_TYPE rows, IDX_TYPE cols, const FP_TYPE& val);
	GlobalDenseMatrix (const GlobalDenseMatrix &src);
	virtual ~GlobalDenseMatrix() {};

    /// return a start index of the active data range
	IDX_TYPE getRangeBegin() const { return static_cast<IDX_TYPE>(0); }

    /// return an end index of the active data range
	IDX_TYPE getRangeEnd() const { return this->_n_rows; }

	/**
	 * Method setZero() set all matrix entries to zero.
	 */
	virtual void setZero();
	/**
	 * Set the value of the matrix entry (row,col) to val.
	 * @param row The index of the row of the matrix.
	 * @param col The index of the column of the matrix.
	 * @param val The value that shoud be set.
	 * @return False if row index or column index are to large, else true.
	 */
	virtual bool setValue(IDX_TYPE row, IDX_TYPE col, FP_TYPE val);

	/**
	 * Method adds a value to the entry at position (row,col).
	 * @param row The index of the row of the matrix.
	 * @param col The index of the column of the matrix.
	 * @param val The value that shoud be added.
	 * @return False if row index or column index are to large, else true.
	 */
	virtual bool add(IDX_TYPE row, IDX_TYPE col, FP_TYPE val);


	/// Add sub-matrix at positions \c row_pos and same column positions as the
	/// given row positions.
	template<class T_DENSE_MATRIX>
	void add(std::vector<IDX_TYPE> const& row_pos,
			const T_DENSE_MATRIX &sub_matrix,
			FP_TYPE fkt = static_cast<FP_TYPE>(1.0))
	{
		this->add(row_pos, row_pos, sub_matrix, fkt);
	}

	template<class T_DENSE_MATRIX>
	void add(RowColumnIndices<IDX_TYPE> const& indices,
			const T_DENSE_MATRIX &sub_matrix,
			FP_TYPE fkt = static_cast<FP_TYPE>(1.0))
	{
		this->add(indices.rows, indices.columns, sub_matrix, fkt);
	}

	template<class T_DENSE_MATRIX>
	void add(std::vector<IDX_TYPE> const& row_pos,
			std::vector<IDX_TYPE> const& col_pos, const T_DENSE_MATRIX &sub_matrix,
			FP_TYPE fkt = static_cast<FP_TYPE>(1.0));

    /// y = mat * x
    void matvec ( const DenseVector<FP_TYPE> &x, DenseVector<FP_TYPE> &y) const;
};

} // end namespace MathLib

#include "GlobalDenseMatrix-impl.h"

#endif /* GLOBALDENSEMATRIX_H_ */
