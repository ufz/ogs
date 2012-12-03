/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file ILinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef ILINEAREQUATION_H_
#define ILINEAREQUATION_H_

#include <iostream>
#include <vector>
#include "BaseLib/Options.h"
#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Sparse/Sparsity.h"

namespace MathLib
{

/**
 * \brief Interface to any algebraic linear equation systems Ax=b
 *
 * This class defines unified interface to algebraic linear equation systems
 * including linear solvers.
 *
 */
class ILinearEquation
{
public:
    /// Index representing invalid
    const static std::size_t index_npos = -1;

    /// 
    virtual ~ILinearEquation()
    {
    }
    ;

    /**
     * create a linear equation 
     *
     * \param dimension    dimension of the equation
     * \param sparsity     sparse pattern
     */
    virtual void create(size_t dimension,
            RowMajorSparsity* sparsity = NULL) = 0;

    /**
     * return if this equation is already created
     */
    virtual bool isCreated() const = 0;

    /**
     * set properties
     *
     * @param option
     */
    virtual void setOption(const BaseLib::Options &option) = 0;

    /**
     * reset the equation
     */
    virtual void reset() = 0;

    /// return dimension of this equation
    virtual size_t getDimension() const = 0;

    /**
     * get an entry in a matrix
     *
     * @param rowId
     * @param colId
     * @return value
     */
    virtual double getA(size_t rowId, size_t colId) const = 0;

    /**
     * set an entry in a matrix
     *
     * @param rowId
     * @param colId
     * @param v
     */
    virtual void setA(size_t rowId, size_t colId, double v) = 0;

    /**
     * add a value into a matrix
     * @param rowId
     * @param colId
     * @param v
     */
    virtual void addA(size_t rowId, size_t colId, double v) = 0;

    /**
     * add a sub-matrix into a matrix
     *
     * @param vec_row_pos   A list of global row index
     * @param vec_col_pos   A list of global column index
     * @param sub_matrix    A sub-matrix
     * @param fkt           scaling parameter
     */
    template <class T_DENSE_MATRIX>
    void addAsub(const std::vector<size_t> &vec_row_pos,
            const std::vector<size_t> &vec_col_pos,
            const T_DENSE_MATRIX &sub_matrix, double fkt = 1.0);

    /**
     * add a sub-matrix into a matrix
     * @param vec_row_col_pos   A list of global row and column index
     * @param sub_matrix        A sub-matrix
     * @param fkt               scaling parameter
     */
    template <class T_DENSE_MATRIX>
    void addAsub(const std::vector<size_t> &vec_row_col_pos,
            const T_DENSE_MATRIX &sub_matrix, double fkt = 1.0);

    /**
     * get RHS entry
     *
     * @param rowId
     * @return
     */
    virtual double getRHS(size_t rowId) const = 0;

    /**
     * get RHS vector
     * @return
     */
    virtual double* getRHS() = 0;

    /**
     * set RHS entry
     *
     * @param rowId
     * @param v
     */
    virtual void setRHS(size_t rowId, double v) = 0;

    /**
     * add RHS entry
     *
     * @param rowId
     * @param v
     */
    virtual void addRHS(size_t rowId, double v) = 0;

    /**
     * add a sub vector to RHS
     *
     * @param vec_row_pos   A list of global row index
     * @param sub_vector    Pointer to a sub-vector
     * @param fkt           Scaling parameter
     */
    inline virtual void addRHSsub(const std::vector<size_t> &vec_row_pos,
            const double* sub_vector, double fkt = 1.0);

    /**
     * add a sub vector to RHS
     *
     * @param vec_row_pos   A list of global row index
     * @param sub_vector    A sub-vector
     * @param fkt           Scaling parameter
     */
    template <class T_DENSE_VECTOR>
    void addRHSsub(const std::vector<size_t> &vec_row_pos,
            const T_DENSE_VECTOR &sub_vector, double fkt = 1.0);

    /**
     * get a solution vector
     *
     * @return a pointer to the vector
     */
    virtual double* getX() = 0;

    /**
     * get an entry in a solution vector
     * @param rowId
     * @return
     */
    virtual double getX(size_t rowId) = 0;

    /**
     * set an entry in a solution vector
     * @param rowId
     * @param v
     */
    virtual void setX(size_t rowId, double v) = 0;

    /**
     * set prescribed value
     * @param row_id
     * @param x
     */
    virtual void setKnownX(size_t row_id, double x) = 0;

    /**
     * set prescribed value
     * @param vec_id    A list of global row index
     * @param vec_x     A list of prescribed values
     */
    virtual void setKnownX(const std::vector<size_t> &vec_id,
            const std::vector<double> &vec_x) = 0;

    /// solve the linear system
    virtual void solve() = 0;

    /// print out the equation for debugging
    virtual void printout(std::ostream &os = std::cout) const = 0;
};

}

#include "ILinearEquation.tpp"

#endif //ILINEAREQUATION_H_
