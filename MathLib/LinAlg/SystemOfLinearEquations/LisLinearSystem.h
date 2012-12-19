/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file LisLinearSystem.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef LISLINEARSYSTEM_H_
#define LISLINEARSYSTEM_H_

#include <iostream>
#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

#include "lis.h"

#include "LisOption.h"

namespace MathLib
{

/**
 * \brief Linear system using Lis solver (http://www.ssisc.org/lis/)
 *
 * This class utilizes Lis solver.
 */
class LisLinearSystem
{
public:
    //---------------------------------------------------------------
    // realization of ISystemOfLinearEquations
    //---------------------------------------------------------------

    /**
     * Create a linear system using Lis
     *
     * @param length    System dimension
     */
    LisLinearSystem(std::size_t length);

    /**
     *
     */
    virtual ~LisLinearSystem();

    /**
     * configure linear solvers
     * @param option
     */
    virtual void setOption(const boost::property_tree::ptree &option);

    /// return the system dimension
    virtual std::size_t getDimension() const { return _dim; };

    /// reset this equation
    virtual void setZero();

    /// set entry in A
    virtual void setMatEntry(std::size_t rowId, std::size_t colId, double v)
    {
        if (rowId==colId)
            _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
        lis_matrix_set_value(LIS_INS_VALUE, rowId, colId, v, _AA);
    }

    /// add value into A
    virtual void addMatEntry(std::size_t rowId, std::size_t colId, double v)
    {
        if (rowId==colId)
            _max_diag_coeff = std::max(_max_diag_coeff, std::abs(v));
        lis_matrix_set_value(LIS_ADD_VALUE, rowId, colId, v, _AA);
    }

    /// get RHS entry
    virtual double getRHSVec(std::size_t rowId) const
    {
        double v;
        lis_vector_get_value(_bb, rowId, &v);
        return v;
    }

    /// set RHS entry
    virtual void setRHSVec(std::size_t rowId, double v)
    {
        lis_vector_set_value(LIS_INS_VALUE, rowId, v, _bb);
    }

    /// add RHS entry
    virtual void addRHSVec(std::size_t rowId, double v)
    {
        lis_vector_set_value(LIS_ADD_VALUE, rowId, v, _bb);
    }

    /// get an entry in a solution vector
    virtual double getSolVec(std::size_t rowId)
    {
        double v;
        lis_vector_get_value(_xx, rowId, &v);
        return v;
    }

    /// set a solution vector
    virtual void setSolVec(std::size_t rowId, double v)
    {
        lis_vector_set_value(LIS_INS_VALUE, rowId, v, _xx);
    }

    /// set prescribed value
    virtual void setKnownSolution(std::size_t rowId, double x)
    {
        _vec_knownX_id.push_back(rowId);
        _vec_knownX_x.push_back(x);
    }

    /// set prescribed values
    virtual void setKnownSolution(const std::vector<std::size_t> &vec_id, const std::vector<double> &vec_x)
    {
        _vec_knownX_id.insert(_vec_knownX_id.end(), vec_id.begin(), vec_id.end());
        _vec_knownX_x.insert(_vec_knownX_x.end(), vec_x.begin(), vec_x.end());
    }

    /// solve this equation
    virtual void solve();

    /// printout this equation for debugging
    virtual void printout(std::ostream &os=std::cout) const;

    //---------------------------------------------------------------
    // specific to this class
    //---------------------------------------------------------------
    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const LisOption &option)
    {
        _option = option;
    }

    /**
     * get linear solver options
     * @return
     */
    LisOption &getOption()
    {
        return _option;
    }

    /**
     * get a solution vector
     * @param x
     */
    void getSolVec(double* x)
    {
        lis_vector_gather(_xx, x);
    }

private:
    bool checkError(int err);
    void applyKnownSolution();

private:
    const std::size_t _dim;
    double _max_diag_coeff;
    LisOption _option;
    LIS_MATRIX _AA;
    LIS_VECTOR _bb;
    LIS_VECTOR _xx;
    std::vector<std::size_t> _vec_knownX_id;
    std::vector<double> _vec_knownX_x;
};


} // MathLib

#endif //LISLINEARSYSTEM_H_

