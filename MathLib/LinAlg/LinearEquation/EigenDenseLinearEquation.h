/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file EigenDenseLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef EIGENDENSELINEAREQUATION_H_
#define EIGENDENSELINEAREQUATION_H_

#include <vector>
#include <Eigen>

#include "ISystemOfLinearEquations.h"


namespace MathLib
{

/**
 * \brief Dense linear equation using Eigen library
 */
class EigenDenseLinearEquation : public ISystemOfLinearEquations
{
public:
    //---------------------------------------------------------------
    // Matrix and vector type
    //---------------------------------------------------------------
    /// Global matrix type
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
    /// Global vector type
    typedef Eigen::VectorXd VectorType;

    //---------------------------------------------------------------
    // Constructor
    //---------------------------------------------------------------
    /**
     *
     */
    EigenDenseLinearEquation(std::size_t length, RowMajorSparsity* sp=NULL);

    /**
     *
     */
    virtual ~EigenDenseLinearEquation() {};

    //---------------------------------------------------------------
    // realization of ISystemOfLinearEquations
    //---------------------------------------------------------------
    virtual void setZero();

    virtual std::size_t getDimension() const { return _A.rows(); }

    virtual double getMatEntry(std::size_t rowId, std::size_t colId) const
    {
        return _A(rowId, colId);
    }

    virtual void setMatEntry(std::size_t rowId, std::size_t colId, double v)
    {
        _A(rowId, colId) = v;
    }

    virtual void addMatEntry(std::size_t rowId, std::size_t colId, double v)
    {
        _A(rowId, colId) += v;
    }

    virtual  double getRHSVec(std::size_t rowId) const
    {
        return _b[rowId];
    }

    virtual void setRHSVec(std::size_t rowId, double v)
    {
        _b[rowId] = v;
    }

    virtual void addRHSVec(std::size_t rowId, double v)
    {
        _b[rowId] += v;
    }

    virtual double getSolVec(std::size_t rowId)
    {
        return _x[rowId];
    }

    virtual void setSolVec(std::size_t rowId, double v)
    {
        _x[rowId] = v;
    }

    virtual void setKnownSolution(std::size_t row_id, double x);

    virtual void setKnownSolution(const std::vector<std::size_t> &vec_id, const std::vector<double> &vec_x);

    virtual void printout(std::ostream &os=std::cout) const;

    virtual void setOption(const BaseLib::Options &/*option*/) {};

    virtual void solve();


    //---------------------------------------------------------------
    // Equation specific functions
    //---------------------------------------------------------------
    /**
     * resize this equation system
     *
     * @param length
     */
    void resize(std::size_t length);

    /**
     * get a LHS matrix
     * @return
     */
    MatrixType* getMatEntry() { return &_A; }

    /**
     * get a solution vector
     * @return
     */
    VectorType* getSolVec() {return &_x;};

    /**
     * get a RHS vector
     * @return
     */
    VectorType* getRHSVec() {return &_b;};

private:
    MatrixType _A;
    VectorType _b;
    VectorType _x;

};

}

#endif //EIGENDENSELINEAREQUATION_H_

