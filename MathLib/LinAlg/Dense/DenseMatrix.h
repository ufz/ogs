/**
 * \file
 * \author Thomas Fischer and Haibing Shao
 * \date   2011-05-24
 * \brief  Definition of the DenseMatrix class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/Vector3.h"

namespace MathLib {

/**
 * Matrix represents a dense matrix for a numeric data type.
 */
template <typename FP_TYPE, typename IDX_TYPE = std::size_t> class DenseMatrix
{
public:
    typedef FP_TYPE FP_T;
    typedef IDX_TYPE IDX_T;

public:
   /// Dense square matrix constructor.
   explicit DenseMatrix (IDX_TYPE rows) : DenseMatrix(rows, rows) {}

   /// Dense rectangular matrix constructor.
   DenseMatrix (IDX_TYPE rows, IDX_TYPE cols);
   DenseMatrix (IDX_TYPE rows, IDX_TYPE cols, const FP_TYPE& val);
   DenseMatrix (const DenseMatrix &src);
   /**
    * Move constructor.
    * @param src The original DenseMatrix object. After applying the
    * move constructor the object src has no rows and columns anymore.
    *
    */
   DenseMatrix (DenseMatrix &&src);

   virtual ~DenseMatrix ();

   /**
    * Get all entries, which are stored in an array.
   */
   const FP_TYPE *getEntries() const
   {
       return _data;
   }

   /**
    * Get all entries, which are stored in an array.
   */
   const FP_TYPE *data() const
   {
       return _data;
   }

   /**
    * Assignment operator, makes a copy of the internal data of the object.
    * @param rhs The DenseMatrix object to the right side of the assignment symbol.
    */
   DenseMatrix& operator=(DenseMatrix const& rhs);

   /**
    * This is the move assignment operator.
    * @param rhs This is the right hand side of a assignment operation.
    * After applying this operation the object src has no rows and columns anymore.
    */
   DenseMatrix& operator=(DenseMatrix && rhs);

   /**
    * Assignment operator
    */
   DenseMatrix& operator=(FP_TYPE const& v);

   /**
    * \f$ y = \alpha \cdot A x + \beta y\f$
    */
   void axpy (FP_TYPE alpha, const FP_TYPE* x, FP_TYPE beta, FP_TYPE* y) const;

   /**
    * DenseMatrix vector multiplication
    */
   FP_TYPE* operator* (FP_TYPE* const& x) const;
   FP_TYPE* operator* (FP_TYPE const* const& x) const;
   template <typename V> V operator* (V const& x) const;
   MathLib::Vector3 operator*(MathLib::Vector3 const& x) const;

   /**
    * DenseMatrix matrix addition.
    */
   DenseMatrix* operator+ (const DenseMatrix& mat) const;
   /**
    * DenseMatrix matrix subtraction
    */
   DenseMatrix* operator- (const DenseMatrix& mat) const;

   /**
    * DenseMatrix matrix multiplication \f$ C = A \cdot B\f$
    * @param mat the matrix \f$ B \f$
    * @return the matrix \f$ C \f$
    */
   DenseMatrix* operator* (const DenseMatrix& mat) const;

   /**
    * matrix transpose
    * @return the transpose of the matrix
    */
   DenseMatrix* transpose() const; // HB & ZC

   /**
    * transpose the matrix in place
    */
   void transposeInPlace();

   DenseMatrix* getSubMatrix (IDX_TYPE b_row, IDX_TYPE b_col, IDX_TYPE e_row, IDX_TYPE e_col) const;

   /**
    * overwrites values of the matrix with the given sub matrix
    * @param b_row the first row
    * @param b_col the first column
    * @param sub_mat the sub matrix
    */
   void setSubMatrix (IDX_TYPE b_row, IDX_TYPE b_col, const DenseMatrix& sub_mat);

   inline FP_TYPE & operator() (IDX_TYPE row, IDX_TYPE col);
   inline FP_TYPE const& operator() (IDX_TYPE row, IDX_TYPE col) const;

   /**
    * writes the matrix entries into the output stream
    * @param out the output stream
    */
   void write (std::ostream& out) const;

    /**
     * get the number of rows
     * @return the number of rows
     */
    IDX_TYPE getNumberOfRows () const { return _n_rows; }
    /**
     * get the number of columns
     * @return the number of columns
     */
    IDX_TYPE getNumberOfColumns () const { return _n_cols; }

    /**
     * get the number of entries in the matrix
     */
    IDX_TYPE size() const { return _n_rows*_n_cols; }

    /**
     * set the identity matrix
     */
    void setIdentity();

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
   //inline IDX_TYPE address(IDX_TYPE i, IDX_TYPE j) const { return j*rows+i; }
   // zero based addressing, C storage layout
   inline IDX_TYPE address(IDX_TYPE i, IDX_TYPE j) const { return i*_n_cols+j; }

   FP_TYPE *_data;
};

/** overload the output operator for class DenseMatrix */
template <typename FP_TYPE, typename IDX_TYPE>
std::ostream& operator<< (std::ostream &os, const DenseMatrix<FP_TYPE, IDX_TYPE> &mat)
{
    mat.write (os);
    return os;
}

} // end namespace MathLib

#include "DenseMatrix-impl.h"
