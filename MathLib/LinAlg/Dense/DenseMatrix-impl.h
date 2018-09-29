/**
 * @file DenseMatrix-impl.h
 * @author Thomas Fischer and Haibing Shao
 * @date Jun 10, 2013
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <cassert>

namespace MathLib
{

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix(IDX_TYPE rows, IDX_TYPE cols) :
        _n_rows(rows), _n_cols(cols), _data(new FP_TYPE[_n_rows * _n_cols])
{}

template<typename FP_TYPE, typename IDX_TYPE>

DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix(IDX_TYPE rows, IDX_TYPE cols,
        FP_TYPE const& initial_value) :
        _n_rows(rows), _n_cols(cols), _data(new FP_TYPE[_n_rows * _n_cols])
{
    const IDX_TYPE n(_n_rows * _n_cols);
    for (IDX_TYPE k(0); k < n; k++)
        _data[k] = initial_value;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix (const DenseMatrix<FP_TYPE, IDX_TYPE>& src) :
        _n_rows(src.getNumberOfRows ()), _n_cols(src.getNumberOfColumns ()), _data (new FP_TYPE[_n_rows * _n_cols])
{
    std::copy(src._data, src._data+_n_rows*_n_cols, _data);
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::DenseMatrix (DenseMatrix<FP_TYPE, IDX_TYPE> &&src) :
       _n_rows(src.getNumberOfRows()), _n_cols(src.getNumberOfColumns())
{
    src._n_rows = 0;
    src._n_cols = 0;
    _data = src._data;
    src._data = nullptr;
}


template <typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>::~DenseMatrix ()
{
   delete [] _data;
}

template <typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator=(DenseMatrix<FP_TYPE, IDX_TYPE> const& rhs)
{
    if (this == &rhs)
        return *this;

    if (_n_rows != rhs.getNumberOfRows() || _n_cols != rhs.getNumberOfColumns()) {
        std::string msg("DenseMatrix::operator=(DenseMatrix const& rhs), Dimension mismatch, ");
        msg += " left hand side: " + std::to_string(_n_rows) + " x "
                + std::to_string(_n_cols);
        msg += " right hand side: " + std::to_string(rhs.getNumberOfRows()) + " x "
                + std::to_string(rhs.getNumberOfColumns());
        throw std::range_error(msg);
        return *this;
    }

    std::copy(rhs._data, rhs._data + _n_rows * _n_cols, _data);

    return *this;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator=(DenseMatrix && rhs)
{
    _n_rows = rhs._n_rows;
    _n_cols = rhs._n_cols;
    _data = rhs._data;

    rhs._n_rows = 0;
    rhs._n_cols = 0;
    rhs._data = nullptr;
    return *this;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator=(FP_TYPE const& v)
{
    std::fill(this->_data, this->_data + this->_n_rows * this->_n_cols, v);
    return *this;
}

template<typename FP_TYPE, typename IDX_TYPE>
void DenseMatrix<FP_TYPE, IDX_TYPE>::axpy(FP_TYPE alpha, const FP_TYPE* x, FP_TYPE beta,
        FP_TYPE* y) const
{
    for (IDX_TYPE i(0); i < _n_rows; i++) {
        y[i] += beta * y[i];
        for (IDX_TYPE j(0); j < _n_cols; j++) {
            y[i] += alpha * _data[address(i, j)] * x[j];
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
    auto* y(new FP_TYPE[_n_rows]);
    for (IDX_TYPE i(0); i < _n_rows; i++) {
        y[i] = 0.0;
        for (IDX_TYPE j(0); j < _n_cols; j++) {
            y[i] += _data[address(i, j)] * x[j];
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
template <typename V>
V DenseMatrix<FP_TYPE, IDX_TYPE>::operator* (V const& x) const
{
    V y(_n_rows);
    for (IDX_TYPE i(0); i < _n_rows; i++) {
        y[i] = 0.0;
        for (IDX_TYPE j(0); j < _n_cols; j++) {
            y[i] += _data[address(i, j)] * x[j];
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
MathLib::Vector3
DenseMatrix<FP_TYPE, IDX_TYPE>::operator*(MathLib::Vector3 const& x) const
{
    assert(_n_rows>2);

    MathLib::Vector3 y;
    for (IDX_TYPE i(0); i < _n_rows; i++) {
        y[i] = 0.0;
        for (IDX_TYPE j(0); j < _n_cols; j++) {
            y[i] += _data[address(i, j)] * x[j];
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::operator+(const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const
{
    // make sure the two matrices have the same dimension.
    if (_n_rows != mat.getNumberOfRows() || _n_cols != mat.getNumberOfColumns())
        throw std::range_error("DenseMatrix::operator+, illegal matrix size!");

    DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE>(_n_rows, _n_cols));
    for (IDX_TYPE i = 0; i < _n_rows; i++) {
        for (IDX_TYPE j = 0; j < _n_cols; j++) {
            (*y)(i, j) = _data[address(i, j)] + mat(i, j);
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::operator-(const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const
{
    // make sure the two matrices have the same dimension.
    if (_n_rows != mat.getNumberOfRows() || _n_cols != mat.getNumberOfColumns())
        throw std::range_error("DenseMatrix::operator-, illegal matrix size!");

    DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE>(_n_rows, _n_cols));
    for (IDX_TYPE i = 0; i < _n_rows; i++) {
        for (IDX_TYPE j = 0; j < _n_cols; j++) {
            (*y)(i, j) = _data[address(i, j)] - mat(i, j);
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::operator*(const DenseMatrix<FP_TYPE, IDX_TYPE>& mat) const
{
    // make sure the two matrices have the same dimension.
    if (_n_cols != mat.getNumberOfRows())
        throw std::range_error(
                "DenseMatrix::operator*, number of rows and cols should be the same!");

    IDX_TYPE y_cols(mat.getNumberOfColumns());
    DenseMatrix<FP_TYPE, IDX_TYPE>* y(
            new DenseMatrix<FP_TYPE, IDX_TYPE>(_n_rows, y_cols, FP_TYPE(0)));

    for (IDX_TYPE i = 0; i < _n_rows; i++) {
        for (IDX_TYPE j = 0; j < y_cols; j++) {
            for (IDX_TYPE k = 0; k < _n_cols; k++)
                (*y)(i, j) += _data[address(i, k)] * mat(k, j);
        }
    }

    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
DenseMatrix<FP_TYPE, IDX_TYPE>*
DenseMatrix<FP_TYPE, IDX_TYPE>::transpose() const
{
    DenseMatrix<FP_TYPE, IDX_TYPE>* y(new DenseMatrix<FP_TYPE, IDX_TYPE>(_n_cols, _n_rows));

    for (IDX_TYPE i = 0; i < _n_rows; i++) {
        for (IDX_TYPE j = 0; j < _n_cols; j++) {
            (*y)(j, i) = _data[address(i, j)];
        }
    }
    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::transposeInPlace()
{
    if (_n_rows==_n_cols) { // square matrix
        for (IDX_TYPE i = 0; i < _n_rows; i++)
            for (IDX_TYPE j = i+1; j < _n_cols; j++)
                std::swap(_data[address(i, j)], _data[address(j, i)]);
    } else { // non-square matrix
        const DenseMatrix<FP_TYPE, IDX_TYPE> org(*this);
        std::swap(_n_rows, _n_cols);
        for (IDX_TYPE i = 0; i < _n_rows; i++) {
            for (IDX_TYPE j = 0; j < _n_cols; j++) {
                _data[address(i, j)] = org(j, i);
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
        throw std::range_error("DenseMatrix::getSubMatrix() illegal sub matrix");
    if (e_row > _n_rows | e_col > _n_cols)
        throw std::range_error("DenseMatrix::getSubMatrix() illegal sub matrix");

    DenseMatrix<FP_TYPE, IDX_TYPE>* y(
            new DenseMatrix<FP_TYPE, IDX_TYPE>(e_row - b_row, e_col - b_col));
    for (IDX_TYPE i = b_row; i < e_row; i++) {
        for (IDX_TYPE j = b_col; j < e_col; j++) {
            (*y)(i - b_row, j - b_col) = _data[address(i, j)];
        }
    }
    return y;
}

template<typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::setSubMatrix(IDX_TYPE b_row, IDX_TYPE b_col,
        const DenseMatrix<FP_TYPE, IDX_TYPE>& sub_mat)
{
    if (b_row + sub_mat.getNumberOfRows() > _n_rows | b_col + sub_mat.getNumberOfColumns() > _n_cols)
        throw std::range_error("DenseMatrix::setSubMatrix() sub matrix to big");

    for (IDX_TYPE i = 0; i < sub_mat.getNumberOfRows(); i++) {
        for (IDX_TYPE j = 0; j < sub_mat.getNumberOfColumns(); j++) {
            _data[address(i + b_row, j + b_col)] = sub_mat(i, j);
        }
    }
}

template<typename FP_TYPE, typename IDX_TYPE>
FP_TYPE&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator() (IDX_TYPE row, IDX_TYPE col)
{
    assert((row < _n_rows) && (col < _n_cols));
    return _data [address(row,col)];
}


template<typename FP_TYPE, typename IDX_TYPE>
FP_TYPE const&
DenseMatrix<FP_TYPE, IDX_TYPE>::operator() (IDX_TYPE row, IDX_TYPE col) const
{
    assert((row < _n_rows) && (col < _n_cols));
    return _data[address(row, col)];
}

template <typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::write (std::ostream &out) const
{
    out << _n_rows << " " << _n_cols << "\n";
    for (IDX_TYPE i = 0; i < _n_rows; i++) {
        for (IDX_TYPE j = 0; j < _n_cols; j++) {
            out << _data[address(i, j)] << "\t";
        }
        out << "\n";
    }
}

template <typename FP_TYPE, typename IDX_TYPE>
void
DenseMatrix<FP_TYPE, IDX_TYPE>::setIdentity()
{
    (*this) = 0.0;
    const IDX_TYPE n_square_rows = std::min(_n_rows, _n_cols);
    for (IDX_TYPE i=0; i<n_square_rows; i++)
        _data[address(i,i)] = 1.0;
}

template <typename FP_TYPE, typename IDX_TYPE>
FP_TYPE
sqrFrobNrm (const DenseMatrix<FP_TYPE, IDX_TYPE> &mat)
{
    FP_TYPE nrm(static_cast<FP_TYPE>(0));
    IDX_TYPE i, j;
    for (j = 0; j < mat.getNumberOfColumns(); j++)
        for (i = 0; i < mat.getNumberOfRows(); i++)
            nrm += mat(i, j) * mat(i, j);

    return nrm;
}

} // end namespace MathLib
