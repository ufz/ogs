/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-24
 * \brief  Definition of the Matrix class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATRIX_TPP
#define MATRIX_TPP

namespace MathLib {

template<class T> Matrix<T>::Matrix (std::size_t rows, std::size_t cols)
      : MatrixBase(rows, cols), _data (new T[_n_rows*_n_cols])
{}

template<class T> Matrix<T>::Matrix (std::size_t rows, std::size_t cols, T const& initial_value)
		: MatrixBase(rows, cols), _data (new T[_n_rows*_n_cols])
{
	const std::size_t n(_n_rows*_n_cols);
	for (std::size_t k(0); k<n; k++)
		_data[k] = initial_value;
}

template<class T> Matrix<T>::Matrix (const Matrix& src) :
	MatrixBase(src.getNRows (), src.getNCols ()), _data (new T[_n_rows * _n_cols])
{
   for (std::size_t i = 0; i < _n_rows; i++)
      for (std::size_t j = 0; j < _n_cols; j++)
         _data[address(i,j)] = src (i, j);
}

template <class T> Matrix<T>::~Matrix ()
{
   delete [] _data;
}

template<class T> void Matrix<T>::axpy ( T alpha, const T* x, T beta, T* y) const
{
   for (std::size_t i(0); i<_n_rows; i++) {
      y[i] += beta * y[i];
      for (std::size_t j(0); j<_n_cols; j++) {
         y[i] += alpha * _data[address(i,j)] * x[j];
      }
   }
}

template<class T> T* Matrix<T>::operator* (const T *x) const
{
	T *y (new T[_n_rows]);
	for (std::size_t i(0); i < _n_rows; i++) {
		y[i] = 0.0;
		for (std::size_t j(0); j < _n_cols; j++) {
			y[i] += _data[address(i, j)] * x[j];
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator+ (const Matrix<T>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (_n_rows != mat.getNRows() || _n_cols != mat.getNCols())
		throw std::range_error("Matrix::operator+, illegal matrix size!");

	Matrix<T>* y(new Matrix<T> (_n_rows, _n_cols));
	for (std::size_t i = 0; i < _n_rows; i++) {
		for (std::size_t j = 0; j < _n_cols; j++) {
			(*y)(i, j) = _data[address(i, j)] + mat(i, j);
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator- (const Matrix<T>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (_n_rows != mat.getNRows() || _n_cols != mat.getNCols())
		throw std::range_error("Matrix::operator-, illegal matrix size!");

	Matrix<T>* y(new Matrix<T> (_n_rows, _n_cols));
	for (std::size_t i = 0; i < _n_rows; i++) {
		for (std::size_t j = 0; j < _n_cols; j++) {
			(*y)(i, j) = _data[address(i, j)] - mat(i, j);
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator* (const Matrix<T>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (_n_cols != mat.getNRows())
		throw std::range_error(
				"Matrix::operator*, number of rows and cols should be the same!");

	std::size_t y_cols(mat.getNCols());
	Matrix<T>* y(new Matrix<T> (_n_rows, y_cols, T(0)));

	for (std::size_t i = 0; i < _n_rows; i++) {
		for (std::size_t j = 0; j < y_cols; j++) {
			for (std::size_t k = 0; k < _n_cols; k++)
				(*y)(i, j) += _data[address(i, k)] * mat(k, j);
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::transpose() const
{
	Matrix<T>* y(new Matrix<T> (_n_cols, _n_rows));

	for (std::size_t i = 0; i < _n_rows; i++) {
		for (std::size_t j = 0; j < _n_cols; j++) {
//			y->_data[y->address(j, i)] = _data[address(i, j)];
			(*y)(j,i) = _data[address(i, j)];
		}
	}
	return y;
}

template<class T> Matrix<T>* Matrix<T>::getSubMatrix(
		std::size_t b_row, std::size_t b_col,
		std::size_t e_row, std::size_t e_col) const throw (std::range_error)
{
	if (b_row >= e_row | b_col >= e_col)
		throw std::range_error ("Matrix::getSubMatrix() illegal sub matrix");
	if (e_row > _n_rows | e_col > _n_cols)
		throw std::range_error ("Matrix::getSubMatrix() illegal sub matrix");

	Matrix<T>* y(new Matrix<T> (e_row-b_row, e_col-b_col));
	for (std::size_t i=b_row; i<e_row; i++) {
		for (std::size_t j=b_col; j<e_col; j++) {
			(*y)(i-b_row, j-b_col) = _data[address(i, j)];
		}
	}
	return y;
}

template<class T> void Matrix<T>::setSubMatrix(
		std::size_t b_row, std::size_t b_col, const Matrix<T>& sub_mat) throw (std::range_error)
{
	if (b_row + sub_mat.getNRows() > _n_rows | b_col + sub_mat.getNCols() > _n_cols)
		throw std::range_error ("Matrix::setSubMatrix() sub matrix to big");

	for (std::size_t i=0; i<sub_mat.getNRows(); i++) {
		for (std::size_t j=0; j<sub_mat.getNCols(); j++) {
			_data[address(i+b_row, j+b_col)] = sub_mat(i,j);
		}
	}
}

template<class T> T& Matrix<T>::operator() (std::size_t row, std::size_t col)
	throw (std::range_error)
{
   if ( (row >= _n_rows) | ( col >= _n_cols) )
	  throw std::range_error ("Matrix: op() const range error");
   return _data [address(row,col)];
}


template<class T> T& Matrix<T>::operator() (std::size_t row, std::size_t col) const
	throw (std::range_error)
{
   if ( (row >= _n_rows) | ( col >= _n_cols) )
      throw std::range_error ("Matrix: op() const range error");
   return _data [address(row,col)];
}

template <class T> void Matrix<T>::write (std::ostream &out) const
{
	for (std::size_t i = 0; i < _n_rows; i++) {
		for (std::size_t j = 0; j < _n_cols; j++) {
			out << _data[address(i, j)] << "\t";
		}
		out << "\n";
	}
}

template <class T> T sqrFrobNrm (const Matrix<T> &mat)
{
	T nrm ((T)(0));
	std::size_t i,j;
	for (j=0; j<mat.getNCols(); j++)
		for (i=0; i<mat.getNRows(); i++)
			nrm += mat(i,j) * mat(i,j);

	return nrm;
}

/** overload the output operator for class Matrix */
template <class T>
std::ostream& operator<< (std::ostream &os, const Matrix<T> &mat)
{
	mat.write (os);
	return os;
}

} // end namespace MathLib

#endif //MATRIX_TPP

