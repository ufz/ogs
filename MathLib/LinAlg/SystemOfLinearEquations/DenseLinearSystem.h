/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-06-25
 * \brief  Definition of the DenseLinearSystem class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DENSELINEARSYSTEM_H_
#define DENSELINEARSYSTEM_H_

#include <iostream>
#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

namespace MathLib
{

class DenseLinearSystem
{
public:
    //---------------------------------------------------------------
    // realization of ISystemOfLinearEquations
    //---------------------------------------------------------------

    /**
     * Create a linear system
     *
     * @param length    System dimension
     */
    DenseLinearSystem(std::size_t length)
    : _mat(length, length), _rhs(length), _x(length)
    {
    }

    /**
     *
     */
    virtual ~DenseLinearSystem();

    /**
     * configure linear solvers
     * @param option
     */
    virtual void setOption(const boost::property_tree::ptree &option);

    /// return the system dimension
    virtual std::size_t getDimension() const { return _rhs.size(); };

    /// reset this equation
    virtual void setZero();

    /// set entry in A
    virtual void setMatEntry(std::size_t rowId, std::size_t colId, double v)
    {
    }

    /// add value into A
    virtual void addMatEntry(std::size_t rowId, std::size_t colId, double v)
    {
    }

    /// get RHS entry
    virtual double getRHSVec(std::size_t rowId) const
    {
        double v=.0;
        return v;
    }

    /// set RHS entry
    virtual void setRHSVec(std::size_t rowId, double v)
    {
    }

    /// add RHS entry
    virtual void addRHSVec(std::size_t rowId, double v)
    {
    }

    /// get an entry in a solution vector
    virtual double getSolVec(std::size_t rowId)
    {
        double v=.0;
        return v;
    }

    /// set a solution vector
    virtual void setSolVec(std::size_t rowId, double v)
    {
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
    typedef MathLib::Matrix<double> MatrixType;
    typedef std::vector<double> VectorType;
    MatrixType& getMat() {return _mat;};
    VectorType& getRHSVec() {return _rhs;};
    VectorType& getSolVec() {return _x;};

private:
    std::vector<std::size_t> _vec_knownX_id;
    std::vector<double> _vec_knownX_x;
    MatrixType _mat;
    VectorType _rhs;
    VectorType _x;
};


} // MathLib

#endif //DENSELINEARSYSTEM_H_

