/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisMatrix class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISMATRIX_H_
#define LISMATRIX_H_

#include <iostream>
#include <cmath>
#include <vector>

#include "lis.h"
#include "LisOption.h"

namespace MathLib
{
class LisVector;

/**
 * \brief Lis matrix wrapper class
 *
 * Lis matrix only supports a regular matrix, i.e. the number of
 * rows should equal to the number of columns.
 */
class LisMatrix
{
public:
    /**
     * constructor
     * @param length
     * @param mat_type default 1 CRS
     */
    LisMatrix(std::size_t length, LisOption::MatrixType mat_type = LisOption::MatrixType::CRS);

    /**
     *
     */
    virtual ~LisMatrix();

    /// return the number of rows
    std::size_t getNRows() const { return _n_rows; };

    /// return the number of columns
    std::size_t getNCols() const { return getNRows(); };

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return _is; }

    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return _ie; }

    /// reset this matrix with keeping its original dimension
    void setZero();

    /// set entry
    int setValue(std::size_t rowId, std::size_t colId, double v)
    {
        if (rowId==colId)
            _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
        lis_matrix_set_value(LIS_INS_VALUE, rowId, colId, v, _AA);
        return 0;
    }

    /// add value
    int addValue(std::size_t rowId, std::size_t colId, double v)
    {
        if (rowId==colId)
            _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
        lis_matrix_set_value(LIS_ADD_VALUE, rowId, colId, v, _AA);
        return 0;
    }

    /// printout this equation for debugging
    void write(const std::string &filename) const;

    /// get a maximum value in diagonal entries
    double getMaxDiagCoeff() const { return _max_diag_coeff; };

    /// return a raw Lis matrix object
    LIS_MATRIX& getRawMatrix() { return _AA; };

    /// y = mat * x
    void matvec ( const LisVector &x, LisVector &y) const;

    ///
	template <class T_DENSE_MATRIX>
	void addSubMatrix(std::vector<std::size_t> const& row_pos,
			std::vector<std::size_t> const& col_pos, const T_DENSE_MATRIX &sub_matrix,
			double fkt = 1.0);

private:
    std::size_t _n_rows;
    double _max_diag_coeff;
    LisOption::MatrixType _mat_type;
    LIS_MATRIX _AA;
    int _is;
    int _ie;

	// friend functions
	friend void finishMatrixAssembly(LisMatrix &mat);
	friend bool isMatrixAssembled(LisMatrix &mat);
};

template<class T_DENSE_MATRIX>
void
LisMatrix::addSubMatrix(std::vector<std::size_t> const& row_pos, std::vector<std::size_t> const& col_pos,
		const T_DENSE_MATRIX &sub_matrix, double fkt)
{
	if (row_pos.size() != sub_matrix.getNRows() || col_pos.size() != sub_matrix.getNCols())
		return;

	const std::size_t n_rows = row_pos.size();
	const std::size_t n_cols = col_pos.size();
	for (std::size_t i = 0; i < n_rows; i++) {
		const std::size_t row = row_pos[i];
		for (std::size_t j = 0; j < n_cols; j++) {
			const std::size_t col = col_pos[j];
			addValue(row, col, fkt * sub_matrix(i, j));
		}
	}
};

/// finish assembly to make this matrix be ready for use
void finishMatrixAssembly(LisMatrix &mat);

/// return if this matrix is already assembled
bool isMatrixAssembled(LisMatrix &mat);

} // MathLib

#endif //LISMATRIX_H_

