/**
 * \file
 * \author Thomas Fischer and Haibing Shao
 * \date Jun 10, 2013
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cassert>

namespace MathLib
{

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix(IDX_TYPE rows, IDX_TYPE cols) :
        n_rows_(rows), n_cols_(cols), data_(new FP_TYPE[n_rows_ * n_cols_])
{}

template<typename FP_TYPE, typename IDX_TYPE>

DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix(IDX_TYPE rows, IDX_TYPE cols,
        FP_TYPE const& initial_value) :
        n_rows_(rows), n_cols_(cols), data_(new FP_TYPE[n_rows_ * n_cols_])
{
    const IDX_TYPE n(n_rows_ * n_cols_);
    for (IDX_TYPE k(0); k < n; k++)
    {
        data_[k] = initial_value;
    }
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix (const DenseMatrix<FP_TYPE, IDX_TYPE>& src) :
        n_rows_(src.getNumberOfRows ()), n_cols_(src.getNumberOfColumns ()), data_ (new FP_TYPE[n_rows_ * n_cols_])
{
    std::copy(src.data_, src.data_+n_rows_*n_cols_, data_);
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix (DenseMatrix<FP_TYPE, IDX_TYPE> &&src) :
       n_rows_(src.getNumberOfRows()), n_cols_(src.getNumberOfColumns())
{
    src.n_rows_ = 0;
    src.n_cols_ = 0;
    data_ = src.data_;
    src.data_ = nullptr;
}


template <typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::~DenseMatrix ()
{
   delete [] data_;
}

template <typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator=(DenseMatrix<FP_TYPE, IDX_TYPE> const& rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    if (n_rows_ != rhs.getNumberOfRows() || n_cols_ != rhs.getNumberOfColumns()) {
        std::string msg("DenseMatrix::operator=(DenseMatrix const& rhs), Dimension mismatch, ");
        msg += " left hand side: " + std::to_string(n_rows_) + " x "
                + std::to_string(n_cols_);
        msg += " right hand side: " + std::to_string(rhs.getNumberOfRows()) + " x "
                + std::to_string(rhs.getNumberOfColumns());
        throw std::range_error(msg);
        return *this;
    }

    std::copy(rhs.data_, rhs.data_ + n_rows_ * n_cols_, data_);

    return *this;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator=(DenseMatrix && rhs)
{
    n_rows_ = rhs.n_rows_;
    n_cols_ = rhs.n_cols_;
    data_ = rhs.data_;

    rhs.n_rows_ = 0;
    rhs.n_cols_ = 0;
    rhs.data_ = nullptr;
    return *this;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator=(FP_TYPE const& v)
{
    std::fill(this->data_, this->data_ + this->n_rows_ * this->n_cols_, v);
    return *this;
}

template<typename FP_TYPE, typename IDX_TYPE>
void DenseMatrix<FP_TYPE, IDX_TYPE>::axpy(FP_TYPE alpha, const FP_TYPE* x, FP_TYPE beta,
        FP_TYPE* y) const
{
    for (IDX_TYPE i(0); i < n_rows_; i++) {
        y[i] += beta * y[i];
        for (IDX_TYPE j(0); j < n_cols_; j++) {
            y[i] += alpha * data_[address(i, j)] * x[j];
        }
    }
}

template<typename FP_TYPE, typename IDX_TYPE>
FP_TYPE* DenseMatrix<FP_TYPE, IDX_TYPE>::operator* (FP_TYPE* const& x) const
{
    return this->operator*(static_cast<FP_TYPE const*>(x));
}

template<typename FP_TYPE, typename IDX_TYPE>
FP_TYPE* DenseMatrix<FP_TYPE, IDX_TYPE>::operator* (FP_TYPE const* const& x) const
{
    auto* y(new FP_TYPE[n_rows_]);
    for (IDX_TYPE i(0); i < n_rows_; i++) {
        y[i] = 0.0;
        for (IDX_TYPE j(0); j < n_cols_; j++) {
            y[i] += data_[address(i, j)] * x[j];
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
template <typename V>
V DenseMatrix<FP_TYPE, IDX_TYPE>::operator* (V const& x) const
{
    V y(n_rows_);
    for (IDX_TYPE i(0); i < n_rows_; i++) {
        y[i] = 0.0;
        for (IDX_TYPE j(0); j < n_cols_; j++) {
            y[i] += data_[address(i, j)] * x[j];
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
MathLib::Vector3
DenseMatrix<FP_TYPE, IDX_TYPE>::operator*(MathLib::Vector3 const& x) const
{
    assert(n_rows_>2);

    MathLib::Vector3 y;
    for (IDX_TYPE i(0); i < n_rows_; i++) {
        y[i] = 0.0;
        for (IDX_TYPE j(0); j < n_cols_; j++) {
            y[i] += data_[address(i, j)] * x[j];
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::operator+(const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const
{
    // make sure the two matrices have the same dimension.
    if (n_rows_ != mat.getNumberOfRows() || n_cols_ != mat.getNumberOfColumns())
    {
        throw std::range_error("DenseMatrix::operator+, illegal matrix size!");
    }

    DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE>(n_rows_, n_cols_));
    for (IDX_TYPE i = 0; i < n_rows_; i++) {
        for (IDX_TYPE j = 0; j < n_cols_; j++) {
            (*y)(i, j) = data_[address(i, j)] + mat(i, j);
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::operator-(const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const
{
    // make sure the two matrices have the same dimension.
    if (n_rows_ != mat.getNumberOfRows() || n_cols_ != mat.getNumberOfColumns())
    {
        throw std::range_error("DenseMatrix::operator-, illegal matrix size!");
    }

    DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE>(n_rows_, n_cols_));
    for (IDX_TYPE i = 0; i < n_rows_; i++) {
        for (IDX_TYPE j = 0; j < n_cols_; j++) {
            (*y)(i, j) = data_[address(i, j)] - mat(i, j);
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::operator*(const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const
{
    // make sure the two matrices have the same dimension.
    if (n_cols_ != mat.getNumberOfRows())
    {
        throw std::range_error(
            "DenseMatrix::operator*, number of rows and cols should be the "
            "same!");
    }

    IDX_TYPE y_cols(mat.getNumberOfColumns());
    DenseMatrix<FP_TYPE, IDX_TYPE>* y(
            new DenseMatrix<FP_TYPE, IDX_TYPE>(n_rows_, y_cols, FP_TYPE(0)));

    for (IDX_TYPE i = 0; i < n_rows_; i++) {
        for (IDX_TYPE j = 0; j < y_cols; j++) {
            for (IDX_TYPE k = 0; k < n_cols_; k++)
            {
                (*y)(i, j) += data_[address(i, k)] * mat(k, j);
            }
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::transpose() const
{
    DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE>(n_cols_, n_rows_));

    for (IDX_TYPE i = 0; i < n_rows_; i++) {
        for (IDX_TYPE j = 0; j < n_cols_; j++) {
            (*y)(j, i) = data_[address(i, j)];
        }
    }
    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::transposeInPlace()
{
    if (n_rows_==n_cols_) { // square matrix
        for (IDX_TYPE i = 0; i < n_rows_; i++)
        {
            for (IDX_TYPE j = i + 1; j < n_cols_; j++)
            {
                std::swap(data_[address(i, j)], data_[address(j, i)]);
            }
        }
    } else { // non-square matrix
        const DenseMatrix<FP_TYPE, IDX_TYPE> org(*this);
        std::swap(n_rows_, n_cols_);
        for (IDX_TYPE i = 0; i < n_rows_; i++) {
            for (IDX_TYPE j = 0; j < n_cols_; j++) {
                data_[address(i, j)] = org(j, i);
            }
        }
    }

}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::getSubMatrix(
        IDX_TYPE b_row, IDX_TYPE b_col,
        IDX_TYPE e_row, IDX_TYPE e_col) const
{
    if (b_row >= e_row | b_col >= e_col)
    {
        throw std::range_error(
            "DenseMatrix::getSubMatrix() illegal sub matrix");
    }
    if (e_row > n_rows_ | e_col > n_cols_)
    {
        throw std::range_error(
            "DenseMatrix::getSubMatrix() illegal sub matrix");
    }

    DenseMatrix<FP_TYPE, IDX_TYPE>* y(
            new DenseMatrix<FP_TYPE, IDX_TYPE>(e_row - b_row, e_col - b_col));
    for (IDX_TYPE i = b_row; i < e_row; i++) {
        for (IDX_TYPE j = b_col; j < e_col; j++) {
            (*y)(i - b_row, j - b_col) = data_[address(i, j)];
        }
    }
    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::setSubMatrix(IDX_TYPE b_row, IDX_TYPE b_col,
        const DenseMatrix<FP_TYPE, IDX_TYPE>& sub_mat)
{
    if (b_row + sub_mat.getNumberOfRows() > n_rows_ |
        b_col + sub_mat.getNumberOfColumns() > n_cols_)
    {
        throw std::range_error("DenseMatrix::setSubMatrix() sub matrix to big");
    }

    for (IDX_TYPE i = 0; i < sub_mat.getNumberOfRows(); i++) {
        for (IDX_TYPE j = 0; j < sub_mat.getNumberOfColumns(); j++) {
            data_[address(i + b_row, j + b_col)] = sub_mat(i, j);
        }
    }
}

template<typename FP_TYPE, typename IDX_TYPE>
FP_TYPE&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator() (IDX_TYPE row, IDX_TYPE col)
{
    assert((row < n_rows_) && (col < n_cols_));
    return data_ [address(row,col)];
}


template<typename FP_TYPE, typename IDX_TYPE>
FP_TYPE const&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator() (IDX_TYPE row, IDX_TYPE col) const
{
    assert((row < n_rows_) && (col < n_cols_));
    return data_[address(row, col)];
}

template <typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::write (std::ostream &out) const
{
    out << n_rows_ << " " << n_cols_ << "\n";
    for (IDX_TYPE i = 0; i < n_rows_; i++) {
        for (IDX_TYPE j = 0; j < n_cols_; j++) {
            out << data_[address(i, j)] << "\t";
        }
        out << "\n";
    }
}

template <typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::setIdentity()
{
    (*this) = 0.0;
    const IDX_TYPE n_square_rows = std::min(n_rows_, n_cols_);
    for (IDX_TYPE i = 0; i < n_square_rows; i++)
    {
        data_[address(i, i)] = 1.0;
    }
}

template <typename FP_TYPE, typename IDX_TYPE>
FP_TYPE
sqrFrobNrm (const DenseMatrix<FP_TYPE, IDX_TYPE> &mat)
{
    FP_TYPE nrm(static_cast<FP_TYPE>(0));
    IDX_TYPE i, j;
    for (j = 0; j < mat.getNumberOfColumns(); j++)
    {
        for (i = 0; i < mat.getNumberOfRows(); i++)
        {
            nrm += mat(i, j) * mat(i, j);
        }
    }

    return nrm;
}

} // end namespace MathLib
