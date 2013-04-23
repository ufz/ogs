/**
 * \file
 * \author Thomas Fischer
 * \date   2011-09-27
 * \brief  Definition of the MatrixBase class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATRIXBASE_H_
#define MATRIXBASE_H_

#include <vector>

namespace MathLib
{

/**
 * class MatrixBase is the basis for all matrix classes (dense and sparse)
 *
 * @tparam FP_TYPE   value type
 * @tparam IDX_TYPE  index type
 */
template <typename FP_TYPE, typename IDX_TYPE>
class MatrixBase
{
public:
	/**
	 * Constructor for initialization of the number of rows and columns
	 * @param nrows number of rows
	 * @param ncols number of columns
	 * @return
	 */
	MatrixBase(unsigned nrows=0, unsigned ncols=0) :
		_n_rows(nrows), _n_cols(ncols)
	{}

	/**
	 * copy constructor.
	 * @param original the object that is copied
	 * @return
	 */
	MatrixBase (MatrixBase<FP_TYPE, IDX_TYPE> const& original) :
		_n_rows (original._n_rows), _n_cols (original._n_cols)
	{}

	/**
	 * destructor of the class.
	 * @return
	 */
	virtual ~MatrixBase() {};
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
     * set all matrix entries to zero
     */
    virtual void setZero() = 0;

    /**
     * set a value to a matrix entry
     *
     * @param rowId
     * @param colId
     * @param v
     * @return 0: if the given matrix entry exists, >0 : else
     */
    virtual int setValue(IDX_TYPE rowId, IDX_TYPE colId, FP_TYPE v) = 0;

    /**
     * add a value to a matrix entry
     *
     * @param rowId
     * @param colId
     * @param v
     * @return 0: if the given matrix entry exists, >0 : else
     */
    virtual int addValue(IDX_TYPE rowId, IDX_TYPE colId, FP_TYPE v) = 0;

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
     * add a sub matrix
     *
     * @param vec_row_pos
     * @param vec_col_pos
     * @param sub_matrix
     * @param fkt
     */
	template <class T_DENSE_MATRIX>
	void addSubMatrix(const std::vector<IDX_TYPE> &vec_row_pos, const std::vector<IDX_TYPE> &vec_col_pos, const T_DENSE_MATRIX &sub_matrix, double fkt=1.0)
	{
	    const std::size_t n_rows = vec_row_pos.size();
	    const std::size_t n_cols = vec_col_pos.size();
	    for (std::size_t i=0; i<n_rows; i++) {
	        const IDX_TYPE rowId = vec_row_pos[i];
	        for (std::size_t j=0; j<n_cols; j++) {
	            const IDX_TYPE colId = vec_col_pos[j];
	            addValue(rowId, colId, fkt*sub_matrix(i,j));
	        }
	    }
	}

protected:
	/**
	 * the number of rows
	 */
	unsigned _n_rows;
	/**
	 * the number of columns
	 */
	unsigned _n_cols;
};

}

#endif /* MATRIXBASE_H_ */
