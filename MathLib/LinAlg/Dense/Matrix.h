/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Matrix.h
 *
 * Created on 2011-05-24 by Thomas Fischer
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <new>
#include <exception>
#include <stdexcept>
#include <iostream>

namespace MathLib {

/**
 * Matrix represents a dense matrix for a numeric data type.
 */
template <class T> class Matrix
{
public:
   Matrix (std::size_t rows, std::size_t cols);
   Matrix (std::size_t rows, std::size_t cols, const T& val);
   Matrix (const Matrix &src);

   ~Matrix ();

   std::size_t getNRows () const { return nrows; }
   std::size_t getNCols () const { return ncols; }
   /**
    * \f$ y = \alpha \cdot A x + \beta y\f$
    */
   void axpy ( T alpha, const T* x, T beta, T* y) const;

   /**
    * Matrix vector multiplication
    * @param x
    * @return
    */
   T* operator* (const T *x) const;
   /**
    * Matrix matrix addition.
    * @param mat
    * @return
    */
   Matrix<T>* operator+ (const Matrix<T>& mat) const throw (std::range_error);
   /**
    * Matrix matrix subtraction
    * @param mat
    * @return
    */
   Matrix<T>* operator- (const Matrix<T>& mat) const throw (std::range_error);

   /**
    * Matrix matrix multiplication \f$ C = A \cdot B\f$
    * @param mat the matrix \f$ B \f$
    * @return the matrix \f$ C \f$
    */
   Matrix<T>* operator* (const Matrix<T>& mat) const throw (std::range_error);

   /**
    * matrix transpose
    * @return the transpose of the matrix
    */
   Matrix<T>* transpose() const; // HB & ZC

   Matrix<T>* getSubMatrix (std::size_t b_row, std::size_t b_col, std::size_t e_row, std::size_t e_col) const throw (std::range_error);

   /**
    * overwrites values of the matrix with the given sub matrix
    * @param b_row the first row
    * @param b_col the first column
    * @param sub_mat the sub matrix
    */
   void setSubMatrix (std::size_t b_row, std::size_t b_col, const Matrix<T>& sub_mat) throw (std::range_error);

   inline T & operator() (std::size_t row, std::size_t col) throw (std::range_error);
   inline T & operator() (std::size_t row, std::size_t col) const throw (std::range_error);

   /**
    * writes the matrix entries into the output stream
    * @param out the output stream
    */
   void write (std::ostream& out) const;

   T const* getCoords () { return data; }

private:
   // zero based addressing, but Fortran storage layout
   //inline std::size_t address(std::size_t i, std::size_t j) const { return j*rows+i; };
   // zero based addressing, C storage layout
   inline std::size_t address(std::size_t i, std::size_t j) const { return i*ncols+j; };

   std::size_t nrows;
   std::size_t ncols;
   T *data;
};

template<class T> Matrix<T>::Matrix (std::size_t rows, std::size_t cols)
      : nrows (rows), ncols (cols), data (new T[nrows*ncols])
{}

template<class T> Matrix<T>::Matrix (std::size_t rows, std::size_t cols, T const& initial_value)
	: nrows (rows), ncols (cols), data (new T[nrows * ncols])
{
	const std::size_t n(nrows*ncols);
	for (std::size_t k(0); k<n; k++) data[k] = initial_value;
}

template<class T> Matrix<T>::Matrix (const Matrix& src) :
	nrows (src.getNRows ()), ncols (src.getNCols ()), data (new T[nrows * ncols])
{
   for (std::size_t i = 0; i < nrows; i++)
      for (std::size_t j = 0; j < ncols; j++)
         data[address(i,j)] = src (i, j);
}

template <class T> Matrix<T>::~Matrix ()
{
   delete [] data;
}

template<class T> void Matrix<T>::axpy ( T alpha, const T* x, T beta, T* y) const
{
   for (std::size_t i(0); i<nrows; i++) {
      y[i] += beta * y[i];
      for (std::size_t j(0); j<ncols; j++) {
         y[i] += alpha * data[address(i,j)] * x[j];
      }
   }
}

template<class T> T* Matrix<T>::operator* (const T *x) const
{
	T *y (new T[nrows]);
	for (std::size_t i(0); i < nrows; i++) {
		y[i] = 0.0;
		for (std::size_t j(0); j < ncols; j++) {
			y[i] += data[address(i, j)] * x[j];
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator+ (const Matrix<T>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (nrows != mat.getNRows() || ncols != mat.getNCols())
		throw std::range_error("Matrix::operator+, illegal matrix size!");

	Matrix<T>* y(new Matrix<T> (nrows, ncols));
	for (std::size_t i = 0; i < nrows; i++) {
		for (std::size_t j = 0; j < ncols; j++) {
			(*y)(i, j) = data[address(i, j)] + mat(i, j);
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator- (const Matrix<T>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (nrows != mat.getNRows() || ncols != mat.getNCols())
		throw std::range_error("Matrix::operator-, illegal matrix size!");

	Matrix<T>* y(new Matrix<T> (nrows, ncols));
	for (std::size_t i = 0; i < nrows; i++) {
		for (std::size_t j = 0; j < ncols; j++) {
			(*y)(i, j) = data[address(i, j)] - mat(i, j);
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::operator* (const Matrix<T>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (ncols != mat.getNRows())
		throw std::range_error(
				"Matrix::operator*, number of rows and cols should be the same!");

	std::size_t y_cols(mat.getNCols());
	Matrix<T>* y(new Matrix<T> (nrows, y_cols, T(0)));

	for (std::size_t i = 0; i < nrows; i++) {
		for (std::size_t j = 0; j < y_cols; j++) {
			for (std::size_t k = 0; k < ncols; k++)
				(*y)(i, j) += data[address(i, k)] * mat(k, j);
		}
	}

	return y;
}

// HS initial implementation
template<class T> Matrix<T>* Matrix<T>::transpose() const
{
	Matrix<T>* y(new Matrix<T> (ncols, nrows));

	for (std::size_t i = 0; i < nrows; i++) {
		for (std::size_t j = 0; j < ncols; j++) {
//			y->data[y->address(j, i)] = data[address(i, j)];
			(*y)(j,i) = data[address(i, j)];
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
	if (e_row > nrows | e_col > ncols)
		throw std::range_error ("Matrix::getSubMatrix() illegal sub matrix");

	Matrix<T>* y(new Matrix<T> (e_row-b_row, e_col-b_col));
	for (std::size_t i=b_row; i<e_row; i++) {
		for (std::size_t j=b_col; j<e_col; j++) {
			(*y)(i-b_row, j-b_col) = data[address(i, j)];
		}
	}
	return y;
}

template<class T> void Matrix<T>::setSubMatrix(
		std::size_t b_row, std::size_t b_col, const Matrix<T>& sub_mat) throw (std::range_error)
{
	if (b_row + sub_mat.getNRows() > nrows | b_col + sub_mat.getNCols() > ncols)
		throw std::range_error ("Matrix::setSubMatrix() sub matrix to big");

	for (std::size_t i=0; i<sub_mat.getNRows(); i++) {
		for (std::size_t j=0; j<sub_mat.getNCols(); j++) {
			data[address(i+b_row, j+b_col)] = sub_mat(i,j);
		}
	}
}

template<class T> T& Matrix<T>::operator() (std::size_t row, std::size_t col)
	throw (std::range_error)
{
   if ( (row >= nrows) | ( col >= ncols) )
	  throw std::range_error ("Matrix: op() const range error");
   return data [address(row,col)];
}


template<class T> T& Matrix<T>::operator() (std::size_t row, std::size_t col) const
	throw (std::range_error)
{
   if ( (row >= nrows) | ( col >= ncols) )
      throw std::range_error ("Matrix: op() const range error");
   return data [address(row,col)];
}

template <class T> void Matrix<T>::write (std::ostream &out) const
{
	for (std::size_t i = 0; i < nrows; i++) {
		for (std::size_t j = 0; j < ncols; j++) {
			out << data[address(i, j)] << "\t";
		}
		out << std::endl;
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

#endif
