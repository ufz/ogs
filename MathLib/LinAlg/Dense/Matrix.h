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

#ifndef MATRIX_H
#define MATRIX_H

#include <new>
#include <exception>
#include <stdexcept>
#include <iostream>

// MathLib/LinAlg
#include "../MatrixBase.h"

namespace MathLib {

/**
 * Matrix represents a dense matrix for a numeric data type.
 */
template <class T> class Matrix : public MatrixBase<T, std::size_t>
{
    using MatrixBase<T, std::size_t>::_n_rows;
    using MatrixBase<T, std::size_t>::_n_cols;
public:
   Matrix (std::size_t rows, std::size_t cols);
   Matrix (std::size_t rows, std::size_t cols, const T& val);
   Matrix (const Matrix &src);

   virtual ~Matrix ();

   /**
    * \f$ y = \alpha \cdot A x + \beta y\f$
    */
   void axpy ( T alpha, const T* x, T beta, T* y) const;

   /**
    * \f$ y = \alpha \cdot A x + \beta y\f$
    */
   void axpy ( T alpha, const Vector<T> &x, T beta, Vector<T> &y) const;

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

   virtual void setZero();

   virtual int setValue(std::size_t row, std::size_t col, T v);

   virtual int addValue(std::size_t row, std::size_t col, T v);

   /**
    * writes the matrix entries into the output stream
    * @param out the output stream
    */
   void write (std::ostream& out) const;

private:
   // zero based addressing, but Fortran storage layout
   //inline std::size_t address(std::size_t i, std::size_t j) const { return j*rows+i; };
   // zero based addressing, C storage layout
   inline std::size_t address(std::size_t i, std::size_t j) const { return i*_n_cols+j; };

   T *_data;
};

} // end namespace MathLib

#include "Matrix.tpp"

#endif
