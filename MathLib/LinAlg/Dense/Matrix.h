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
#include <vector>

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

   virtual ~Matrix ();

   /**
    * get the number of rows
    * @return the number of rows
    */
   unsigned getNRows () const { return _n_rows; }

   /**
    * get the number of columns
    * @return the number of columns
    */
   unsigned getNCols () const { return _n_cols; }

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

   /**
    * add values of a sub matrix into the matrix
    *
    * @param vec_row_pos    positions of each row of sub-matrix in this matrix
    * @param vec_col_pos    positions of each column of sub-matrix entry in this matrix
    * @param sub_matrix     sub-matrix
    * @param factor
    */
   template <class T_DENSE_MATRIX>
   void addSubMatrix(const std::vector<std::size_t> &vec_row_pos, const std::vector<std::size_t> &vec_col_pos, const T_DENSE_MATRIX &sub_matrix, double factor=1.0)
   {
       const std::size_t n_rows = vec_row_pos.size();
       const std::size_t n_cols = vec_col_pos.size();
       for (std::size_t i=0; i<n_rows; i++) {
           const std::size_t rowId = vec_row_pos[i];
           for (std::size_t j=0; j<n_cols; j++) {
               const std::size_t colId = vec_col_pos[j];
               addValue(rowId, colId, factor*sub_matrix(i,j));
           }
       }
   }

   inline T & operator() (std::size_t row, std::size_t col) throw (std::range_error);
   inline T & operator() (std::size_t row, std::size_t col) const throw (std::range_error);

   virtual void setZero();

   virtual int setValue(std::size_t row, std::size_t col, T v);

   virtual int addValue(std::size_t row, std::size_t col, T v);

   /**
    * complete the matrix assembly
    *
    * This function should be called before using the matrix because some matrix
    * implementations store results of setValue() and addValue() in caches.
    */
   virtual void finishAssembly() {};

   /**
    * return if this matrix is ready to use
    *
    * @return
    */
   virtual bool isAssembled() const { return true; };

   /**
    * writes the matrix entries into the output stream
    * @param out the output stream
    */
   void write (std::ostream& out) const;

private:
   // zero based addressing, but Fortran storage layout
   //inline std::size_t address(std::size_t i, std::size_t j) const { return j*rows+i; };
   // zero based addressing, C storage layout
   inline std::size_t address(std::size_t i, std::size_t j) const { return i*this->_n_cols+j; };

   unsigned _n_rows;
   unsigned _n_cols;
   T *_data;
};

} // end namespace MathLib

#include "Matrix.tpp"

#endif
