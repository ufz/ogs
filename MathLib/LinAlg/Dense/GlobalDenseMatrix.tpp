/**
 * @file GlobalDenseMatrix.tpp
 * @author Thomas Fischer
 * @date Jun 10, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

namespace MathLib
{

template<typename FP_TYPE, typename IDX_TYPE>
GlobalDenseMatrix<FP_TYPE, IDX_TYPE>::GlobalDenseMatrix(IDX_TYPE rows, IDX_TYPE cols) :
		DenseMatrix<FP_TYPE, IDX_TYPE>(rows, cols)
{}

template<typename FP_TYPE, typename IDX_TYPE>
GlobalDenseMatrix<FP_TYPE, IDX_TYPE>::GlobalDenseMatrix (IDX_TYPE rows, IDX_TYPE cols, const FP_TYPE& val) :
		DenseMatrix<FP_TYPE, IDX_TYPE>(rows, cols, val)
{}

template<typename FP_TYPE, typename IDX_TYPE>
GlobalDenseMatrix<FP_TYPE, IDX_TYPE>::GlobalDenseMatrix(const GlobalDenseMatrix &src) :
		DenseMatrix<FP_TYPE, IDX_TYPE>(src.getNRows(), src.getNCols())
{
	std::copy(src._data, src._data + this->_n_rows * this->_n_cols, this->_data);
}

template<typename FP_TYPE, typename IDX_TYPE>
void
GlobalDenseMatrix<FP_TYPE, IDX_TYPE>::setZero()
	{
		std::fill(this->_data, this->_data+this->_n_rows*this->_n_cols, static_cast<FP_TYPE>(0));
	}

template<typename FP_TYPE, typename IDX_TYPE>
bool
GlobalDenseMatrix<FP_TYPE, IDX_TYPE>::setValue(IDX_TYPE row, IDX_TYPE col, FP_TYPE val)
	{
		if (row >= this->_n_rows || col >= this->_n_cols)
			return false;
		this->operator()(row,col) = val;
		return true;
	}

template<typename FP_TYPE, typename IDX_TYPE>
bool
GlobalDenseMatrix<FP_TYPE, IDX_TYPE>::addValue(IDX_TYPE row, IDX_TYPE col, FP_TYPE val)
{
	if (row >= this->_n_rows || col >= this->_n_cols)
		return false;
	this->operator()(row, col) += val;
	return true;
}

template<typename FP_TYPE, typename IDX_TYPE>
template<class T_DENSE_MATRIX>
void
GlobalDenseMatrix<FP_TYPE, IDX_TYPE>::addSubMatrix(std::vector<IDX_TYPE> const& row_pos, std::vector<IDX_TYPE> const& col_pos,
		const T_DENSE_MATRIX &sub_matrix, FP_TYPE fkt)
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


} // end namespace MathLib


