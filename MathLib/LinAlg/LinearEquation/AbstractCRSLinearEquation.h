/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file AbstractCRSLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef ABSTRACTCRSLINEAREQUATION_H_
#define ABSTRACTCRSLINEAREQUATION_H_

#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/Sparse/CRSMatrix.h"
#include "ISystemOfLinearEquations.h"

namespace MathLib
{

/**
 * \brief Abstract class of any linear equations based on CRS sparse matrix
 *
 * \tparam IDX_TYPE     Index data type, i.e. signed or unsigned
 */
template<typename IDX_TYPE>
class AbstractCRSLinearEquation : public ISystemOfLinearEquations
{
public:
    //---------------------------------------------------------------
    // Matrix and vector type
    //---------------------------------------------------------------
    typedef CRSMatrix<double, IDX_TYPE> MatrixType;
    typedef std::vector<double> VectorType;
    
    //---------------------------------------------------------------
    // Constructor
    //---------------------------------------------------------------
    /**
     * create a linear equation
     *
     * \param dim       dimension of the system
     * \param sp        sparse pattern
     */
    AbstractCRSLinearEquation(size_t length, RowMajorSparsity* sp);

    /// 
    virtual ~AbstractCRSLinearEquation()
    {
        if (_A!=NULL) {
            delete _A;
            _A = NULL;
        }
    }

    //---------------------------------------------------------------
    // realization of ISystemOfLinearEquations
    //---------------------------------------------------------------
    /// return the dimension
    virtual size_t getDimension() const { return _x.size(); };

    /// reset this equation
    virtual void setZero();

    /// get entry in A
    virtual double getMatEntry(size_t rowId, size_t colId) const
    {
        return _A->getValue(rowId, colId);
    }

    /// set entry in A
    virtual void setMatEntry(size_t rowId, size_t colId, double v)
    {
        _A->setValue(rowId, colId, v);
    }

    /// add value into A
    virtual void addMatEntry(size_t rowId, size_t colId, double v)
    {
        _A->addValue(rowId, colId, v);
    }

    /// get RHS entry
    virtual double getRHSVec(size_t rowId) const { return _b[rowId]; }

    /// set RHS entry
    virtual void setRHSVec(size_t rowId, double v) { _b[rowId] = v; }

    /// add RHS entry
    virtual void addRHSVec(size_t rowId, double v)
    {
        _b[rowId] += v;
    }

    /// get an entry in a solution vector
    virtual double getSolVec(size_t rowId)
    {
        return _x[rowId];
    }

    /// set a solution vector
    virtual void setSolVec(size_t i, double v)
    {
        _x[i] = v;
    }

    /// set prescribed value
    virtual void setKnownSolution(size_t id, double x)
    {
        _vec_knownX_id.push_back(id);
        _vec_knownX_x.push_back(x);
    }

    /// set prescribed values
    virtual void setKnownSolution(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
    {
        _vec_knownX_id.insert(_vec_knownX_id.end(), vec_id.begin(), vec_id.end());
        _vec_knownX_x.insert(_vec_knownX_x.end(), vec_x.begin(), vec_x.end());
    }

    /// solve this equation
    virtual void solve();

    /// printout this equation for debugging
    virtual void printout(std::ostream &os=std::cout) const
    {
        os << "#A=" << std::endl;
        _A->printMat();
        os << "#x=" << std::endl;
        for (size_t i=0; i<_x.size(); i++) os << _x[i] << " "; os << std::endl;
        os << "#b=" << std::endl;
        for (size_t i=0; i<_b.size(); i++) os << _b[i] << " "; os << std::endl;
    }

    //---------------------------------------------------------------
    // Equation specific functions
    //---------------------------------------------------------------
    /// return A matrix object
    MatrixType* getMatEntry() { return _A;};

    /// get RHS vector
    double* getRHSVec() { return &_b[0]; }

    /// get a solution vector
    double* getSolVec() { return &_x[0]; }
    
protected:
    virtual void solveEqs(MatrixType *A, double *rhs, double *x) = 0;

private:
    void setKnownXi_ReduceSizeOfEQS(MatrixType *A, double *org_eqsRHS, double *org_eqsX, const std::vector<size_t> &vec_id, const std::vector<double> &vec_x, std::vector<double> &out_b, std::vector<double> &out_x, std::map<size_t,size_t> &map_solved_orgEqs);

    DISALLOW_COPY_AND_ASSIGN(AbstractCRSLinearEquation);

private:
    MatrixType* _A;
    VectorType _b;
    VectorType _x;
    std::vector<size_t> _vec_knownX_id;
    std::vector<double> _vec_knownX_x;
};

} //end

#include "AbstractCRSLinearEquation.tpp"

#endif //ABSTRACTCRSLINEAREQUATION_H_
