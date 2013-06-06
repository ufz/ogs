/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-24
 * \brief  Definition of the DenseMatrix class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <new>
#include <exception>
#include <stdexcept>
#include <iostream>

// BaseLib
#include "StringTools.h"

namespace MathLib {

/**
 * Matrix represents a dense matrix for a numeric data type.
 */
template <typename FP_TYPE, typename IDX_TYPE = std::size_t> class DenseMatrix
{
public:
   DenseMatrix (IDX_TYPE rows, IDX_TYPE cols);
   DenseMatrix (IDX_TYPE rows, IDX_TYPE cols, const FP_TYPE& val);
   DenseMatrix (const DenseMatrix &src);

   virtual ~DenseMatrix ();

   DenseMatrix& operator=(DenseMatrix const& rhs) throw (std::range_error)
   {
	   if (this == &rhs)
		   return *this;

	   if (_n_rows != rhs.getNRows() || _n_cols != rhs.getNCols()) {
		   std::string msg("DenseMatrix::operator=, Dimension mismatch, ");
		   msg += " left hand side: "  + BaseLib::number2str(_n_rows) + " x " + BaseLib::number2str(_n_cols);
		   msg += " right hand side: "  + BaseLib::number2str(rhs.getNRows()) + " x " + BaseLib::number2str(rhs.getNCols());
		   throw std::range_error(msg);
		   return *this;
	   }

	   std::copy(rhs._data, rhs._data+_n_rows*_n_cols, _data);

	   return *this;
   }
   /**
    * \f$ y = \alpha \cdot A x + \beta y\f$
    */
   void axpy (FP_TYPE alpha, const FP_TYPE* x, FP_TYPE beta, FP_TYPE* y) const;

   /**
    * DenseMatrix vector multiplication
    * @param x
    * @return
    */
   FP_TYPE* operator* (const FP_TYPE *x) const;
   /**
    * DenseMatrix matrix addition.
    * @param mat
    * @return
    */
   DenseMatrix<FP_TYPE, IDX_TYPE>* operator+ (const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const throw (std::range_error);
   /**
    * DenseMatrix matrix subtraction
    * @param mat
    * @return
    */
   DenseMatrix<FP_TYPE, IDX_TYPE>* operator- (const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const throw (std::range_error);

   /**
    * DenseMatrix matrix multiplication \f$ C = A \cdot B\f$
    * @param mat the matrix \f$ B \f$
    * @return the matrix \f$ C \f$
    */
   DenseMatrix<FP_TYPE, IDX_TYPE>* operator* (const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const throw (std::range_error);

   /**
    * matrix transpose
    * @return the transpose of the matrix
    */
   DenseMatrix<FP_TYPE, IDX_TYPE>* transpose() const; // HB & ZC

   DenseMatrix<FP_TYPE, IDX_TYPE>* getSubMatrix (IDX_TYPE b_row, IDX_TYPE b_col, IDX_TYPE e_row, IDX_TYPE e_col) const throw (std::range_error);

   /**
    * overwrites values of the matrix with the given sub matrix
    * @param b_row the first row
    * @param b_col the first column
    * @param sub_mat the sub matrix
    */
   void setSubMatrix (IDX_TYPE b_row, IDX_TYPE b_col, const DenseMatrix<FP_TYPE, IDX_TYPE>& sub_mat) throw (std::range_error);

   inline FP_TYPE & operator() (IDX_TYPE row, IDX_TYPE col) throw (std::range_error);
   inline FP_TYPE const& operator() (IDX_TYPE row, IDX_TYPE col) const throw (std::range_error);

   /**
    * writes the matrix entries into the output stream
    * @param out the output stream
    */
   void write (std::ostream& out) const;

	/**
	 * get the number of rows
	 * @return the number of rows
	 */
	IDX_TYPE getNRows () const { return _n_rows; }
	/**
	 * get the number of columns
	 * @return the number of columns
	 */
	IDX_TYPE getNCols () const { return _n_cols; }

protected:
	/**
	 * the number of rows
	 */
	IDX_TYPE _n_rows;
	/**
	 * the number of columns
	 */
	IDX_TYPE _n_cols;

   // zero based addressing, but Fortran storage layout
   //inline IDX_TYPE address(IDX_TYPE i, IDX_TYPE j) const { return j*rows+i; };
   // zero based addressing, C storage layout
   inline IDX_TYPE address(IDX_TYPE i, IDX_TYPE j) const { return i*_n_cols+j; };

   FP_TYPE *_data;
};

template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix (IDX_TYPE rows, IDX_TYPE cols)
      : _n_rows(rows), _n_cols(cols), _data (new FP_TYPE[_n_rows*_n_cols])
{}

template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix (IDX_TYPE rows, IDX_TYPE cols, FP_TYPE const& initial_value)
		: _n_rows(rows), _n_cols(cols), _data (new FP_TYPE[_n_rows*_n_cols])
{
	const IDX_TYPE n(_n_rows*_n_cols);
	for (IDX_TYPE k(0); k<n; k++)
		_data[k] = initial_value;
}

template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix (const DenseMatrix& src) :
		_n_rows(src.getNRows ()), _n_cols(src.getNCols ()), _data (new FP_TYPE[_n_rows * _n_cols])
{
	std::copy(src._data, src._data+_n_rows*_n_cols, _data);
}

template <typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>::~DenseMatrix ()
{
   delete [] _data;
}

template<typename FP_TYPE, typename IDX_TYPE> void DenseMatrix<FP_TYPE, IDX_TYPE>::axpy ( FP_TYPE alpha, const FP_TYPE* x, FP_TYPE beta, FP_TYPE* y) const
{
   for (IDX_TYPE i(0); i<_n_rows; i++) {
      y[i] += beta * y[i];
      for (IDX_TYPE j(0); j<_n_cols; j++) {
         y[i] += alpha * _data[address(i,j)] * x[j];
      }
   }
}

template<typename FP_TYPE, typename IDX_TYPE> FP_TYPE* DenseMatrix<FP_TYPE, IDX_TYPE>::operator* (const FP_TYPE *x) const
{
	FP_TYPE *y (new FP_TYPE[_n_rows]);
	for (IDX_TYPE i(0); i < _n_rows; i++) {
		y[i] = 0.0;
		for (IDX_TYPE j(0); j < _n_cols; j++) {
			y[i] += _data[address(i, j)] * x[j];
		}
	}

	return y;
}

// HS initial implementation
template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>* DenseMatrix<FP_TYPE, IDX_TYPE>::operator+ (const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (_n_rows != mat.getNRows() || _n_cols != mat.getNCols())
		throw std::range_error("DenseMatrix::operator+, illegal matrix size!");

	DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE> (_n_rows, _n_cols));
	for (IDX_TYPE i = 0; i < _n_rows; i++) {
		for (IDX_TYPE j = 0; j < _n_cols; j++) {
			(*y)(i, j) = _data[address(i, j)] + mat(i, j);
		}
	}

	return y;
}

// HS initial implementation
template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>* DenseMatrix<FP_TYPE, IDX_TYPE>::operator- (const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (_n_rows != mat.getNRows() || _n_cols != mat.getNCols())
		throw std::range_error("DenseMatrix::operator-, illegal matrix size!");

	DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE> (_n_rows, _n_cols));
	for (IDX_TYPE i = 0; i < _n_rows; i++) {
		for (IDX_TYPE j = 0; j < _n_cols; j++) {
			(*y)(i, j) = _data[address(i, j)] - mat(i, j);
		}
	}

	return y;
}

// HS initial implementation
template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>* DenseMatrix<FP_TYPE, IDX_TYPE>::operator* (const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const throw (std::range_error)
{
	// make sure the two matrices have the same dimension.
	if (_n_cols != mat.getNRows())
		throw std::range_error(
				"DenseMatrix::operator*, number of rows and cols should be the same!");

	IDX_TYPE y_cols(mat.getNCols());
	DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE> (_n_rows, y_cols, FP_TYPE(0)));

	for (IDX_TYPE i = 0; i < _n_rows; i++) {
		for (IDX_TYPE j = 0; j < y_cols; j++) {
			for (IDX_TYPE k = 0; k < _n_cols; k++)
				(*y)(i, j) += _data[address(i, k)] * mat(k, j);
		}
	}

	return y;
}

// HS initial implementation
template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>* DenseMatrix<FP_TYPE, IDX_TYPE>::transpose() const
{
	DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE> (_n_cols, _n_rows));

	for (IDX_TYPE i = 0; i < _n_rows; i++) {
		for (IDX_TYPE j = 0; j < _n_cols; j++) {
//			y->_data[y->address(j, i)] = _data[address(i, j)];
			(*y)(j,i) = _data[address(i, j)];
		}
	}
	return y;
}

template<typename FP_TYPE, typename IDX_TYPE> DenseMatrix<FP_TYPE, IDX_TYPE>* DenseMatrix<FP_TYPE, IDX_TYPE>::getSubMatrix(
		IDX_TYPE b_row, IDX_TYPE b_col,
		IDX_TYPE e_row, IDX_TYPE e_col) const throw (std::range_error)
{
	if (b_row >= e_row | b_col >= e_col)
		throw std::range_error ("DenseMatrix::getSubMatrix() illegal sub matrix");
	if (e_row > _n_rows | e_col > _n_cols)
		throw std::range_error ("DenseMatrix::getSubMatrix() illegal sub matrix");

	DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE> (e_row-b_row, e_col-b_col));
	for (IDX_TYPE i=b_row; i<e_row; i++) {
		for (IDX_TYPE j=b_col; j<e_col; j++) {
			(*y)(i-b_row, j-b_col) = _data[address(i, j)];
		}
	}
	return y;
}

template<typename FP_TYPE, typename IDX_TYPE> void DenseMatrix<FP_TYPE, IDX_TYPE>::setSubMatrix(
		IDX_TYPE b_row, IDX_TYPE b_col, const DenseMatrix<FP_TYPE, IDX_TYPE>& sub_mat) throw (std::range_error)
{
	if (b_row + sub_mat.getNRows() > _n_rows | b_col + sub_mat.getNCols() > _n_cols)
		throw std::range_error ("DenseMatrix::setSubMatrix() sub matrix to big");

	for (IDX_TYPE i=0; i<sub_mat.getNRows(); i++) {
		for (IDX_TYPE j=0; j<sub_mat.getNCols(); j++) {
			_data[address(i+b_row, j+b_col)] = sub_mat(i,j);
		}
	}
}

template<typename FP_TYPE, typename IDX_TYPE> FP_TYPE& DenseMatrix<FP_TYPE, IDX_TYPE>::operator() (IDX_TYPE row, IDX_TYPE col)
	throw (std::range_error)
{
   if ( (row >= _n_rows) | ( col >= _n_cols) )
	  throw std::range_error ("DenseMatrix: op() const range error");
   return _data [address(row,col)];
}


template<typename FP_TYPE, typename IDX_TYPE> FP_TYPE const& DenseMatrix<FP_TYPE, IDX_TYPE>::operator() (IDX_TYPE row, IDX_TYPE col) const
	throw (std::range_error)
{
   if ( (row >= _n_rows) | ( col >= _n_cols) )
      throw std::range_error ("DenseMatrix: op() const range error");
   return _data [address(row,col)];
}

template <typename FP_TYPE, typename IDX_TYPE> void DenseMatrix<FP_TYPE, IDX_TYPE>::write (std::ostream &out) const
{
	for (IDX_TYPE i = 0; i < _n_rows; i++) {
		for (IDX_TYPE j = 0; j < _n_cols; j++) {
			out << _data[address(i, j)] << "\t";
		}
		out << "\n";
	}
}

template <typename FP_TYPE, typename IDX_TYPE> FP_TYPE sqrFrobNrm (const DenseMatrix<FP_TYPE, IDX_TYPE> &mat)
{
	FP_TYPE nrm ((FP_TYPE)(0));
	IDX_TYPE i,j;
	for (j=0; j<mat.getNCols(); j++)
		for (i=0; i<mat.getNRows(); i++)
			nrm += mat(i,j) * mat(i,j);

	return nrm;
}

/** overload the output operator for class DenseMatrix */
template <typename FP_TYPE, typename IDX_TYPE>
std::ostream& operator<< (std::ostream &os, const DenseMatrix<FP_TYPE, IDX_TYPE> &mat)
{
	mat.write (os);
	return os;
}

} // end namespace MathLib

#endif
