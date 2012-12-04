/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file ISystemOfLinearEquations.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef ISYSTEMOFLINEAREQUATIONS_H_
#define ISYSTEMOFLINEAREQUATIONS_H_

#include <iostream>
#include <vector>
#include "BaseLib/Options.h"
#include "MathLib/LinAlg/Sparse/Sparsity.h"

namespace MathLib
{

/**
 * \brief Interface to any system of linear equations Ax=b
 *
 * This class defines unified interface to a system of linear equations
 * including linear solvers.
 *
 */
class ISystemOfLinearEquations
{
public:
    /// Index representing invalid
    const static std::size_t index_npos = -1;

    /**
     * Constructor
     *
     * \param dimension    system dimension
     * \param sparsity     sparse pattern
     */
    ISystemOfLinearEquations(std::size_t /*dimension*/,
            RowMajorSparsity* /*sparsity*/ = NULL) {};

    /// 
    virtual ~ISystemOfLinearEquations()
    {
    }

    /**
     * set properties
     *
     * @param option
     */
    virtual void setOption(const BaseLib::Options &option) = 0;

    /**
     * initialize this system with zero
     *
     * This function initialize the system by setting zero to all entries in
     * a matrix, a RHS vector and solution vector. Dimension of the system is
     * not changed.
     */
    virtual void setZero() = 0;

    /// return dimension of this equation
    virtual std::size_t getDimension() const = 0;

    /**
     * get an entry in a matrix
     *
     * @param rowId
     * @param colId
     * @return value
     */
    virtual double getMatEntry(std::size_t rowId, std::size_t colId) const = 0;

    /**
     * set an entry in a matrix
     *
     * @param rowId
     * @param colId
     * @param v
     */
    virtual void setMatEntry(std::size_t rowId, std::size_t colId, double v) = 0;

    /**
     * add a value into a matrix
     *
     * @param rowId
     * @param colId
     * @param v
     */
    virtual void addMatEntry(std::size_t rowId, std::size_t colId, double v) = 0;

    /**
     * add a sub-matrix into a matrix at given row and column positions
     *
     * @param vec_row_pos
     *  A vector of global row index. The number of entries of vec_row_pos have
     *  to be identical to the number of rows of the local matrix.
     * @param vec_col_pos
     *  A vector of global column index. The number of vec_col_pos have to be
     *  identical to the number of columns of the local matrix.
     * @param sub_matrix    A sub-matrix
     * @param fac           scaling parameter
     */
    template <class T_DENSE_MATRIX>
    void addSubMat(const std::vector<std::size_t> &vec_row_pos,
            const std::vector<std::size_t> &vec_col_pos,
            const T_DENSE_MATRIX &sub_matrix, double fac = 1.0);

    /**
     * add a sub-matrix into a matrix
     *
     * @param vec_row_col_pos
     *  A vector of global row and column index. The number of entries of
     *  vec_row_col_pos have to be identical to the number of rows and columns
     *  of the local matrix.
     * @param sub_matrix        A sub-matrix
     * @param fac               scaling parameter
     */
    template <class T_DENSE_MATRIX>
    void addSubMat(const std::vector<std::size_t> &vec_row_col_pos,
            const T_DENSE_MATRIX &sub_matrix, double fac = 1.0);

    /**
     * get RHS entry
     *
     * @param rowId
     * @return
     */
    virtual double getRHSVec(std::size_t rowId) const = 0;

    /**
     * set RHS entry
     *
     * @param rowId
     * @param v
     */
    virtual void setRHSVec(std::size_t rowId, double v) = 0;

    /**
     * add RHS entry
     *
     * @param rowId
     * @param v
     */
    virtual void addRHSVec(std::size_t rowId, double v) = 0;

    /**
     * add a sub vector to RHS
     *
     * @param vec_row_pos
     *  A vector of global row index. The number of entries of vec_row_pos
     *  have to be identical to the number of rows of the sub vector.
     * @param sub_vector    Pointer to a sub-vector
     * @param fac           Scaling parameter
     */
    inline virtual void addSubRHS(const std::vector<std::size_t> &vec_row_pos,
            const double* sub_vector, double fac = 1.0);

    /**
     * add a sub vector to RHS
     *
     * @param vec_row_pos
     *  A vector of global row index. The number of entries of vec_row_pos
     *  have to be identical to the number of rows of the sub vector.
     * @param sub_vector    A sub-vector
     * @param fac           Scaling parameter
     */
    template <class T_DENSE_VECTOR>
    void addSubRHS(const std::vector<std::size_t> &vec_row_pos,
            const T_DENSE_VECTOR &sub_vector, double fac = 1.0);

    /**
     * get an entry in a solution vector
     * @param rowId
     * @return
     */
    virtual double getSolVec(std::size_t rowId) = 0;

    /**
     * set an entry in a solution vector
     * @param rowId
     * @param v
     */
    virtual void setSolVec(std::size_t rowId, double v) = 0;

    /**
     * set prescribed value
     * @param row_id
     * @param x
     */
    virtual void setKnownSolution(std::size_t row_id, double x) = 0;

    /**
     * set prescribed value
     * @param vec_id    A vector of global row index
     * @param vec_x     A vector of prescribed values
     */
    virtual void setKnownSolution(const std::vector<std::size_t> &vec_id,
            const std::vector<double> &vec_x) = 0;

    /// solve the linear system
    virtual void solve() = 0;

    /// print out the equation for debugging
    virtual void printout(std::ostream &os = std::cout) const = 0;

};

}

#include "ISystemOfLinearEquations.tpp"

#endif //ISYSTEMOFLINEAREQUATIONS_H_
