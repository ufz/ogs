/**
 * \file   IDiscreteLinearSystem.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef IDISCRETELINEARSYSTEM_H_
#define IDISCRETELINEARSYSTEM_H_

#include <vector>

#include "DiscreteObjectWithID.h"
#include "IDiscreteVector.h"

namespace DiscreteLib
{

// forward declaration
class DofEquationIdTable;
class IDiscreteLinearSystemAssembler;

/** 
 * \brief Interface of discrete linear systems
 */
class IDiscreteLinearSystem: public DiscreteObjectWithID
{
public:
    /// RHS and solution vector type
    typedef IDiscreteVector<double> GlobalVectorType;

    ///
    virtual ~IDiscreteLinearSystem(){};

    /**
     * zero the linear system
     */
    virtual void setZero() = 0;

    /**
     * construct this linear system
     *
     * @param assembler
     */
    virtual void construct(IDiscreteLinearSystemAssembler& assembler) = 0;

    /**
     * set
     * @param dofId
     * @param list_discrete_pt_id
     * @param list_prescribed_values
     */
    virtual void setKnownSolution(
                    std::size_t dofId,
                    const std::vector<std::size_t> &list_discrete_pt_id,
                    const std::vector<double> &list_prescribed_values) = 0;


    /**
     * solve this linear system
     */
    virtual void solve() = 0;

    /**
     * get a DoF mapping table
     * @return
     */
    virtual DofEquationIdTable* getDofEquationIdTable() const = 0;

    /**
     * get a solution vector
     * @param x
     */
    virtual void getSolVec(std::vector<double> &x) = 0;

    /**
     * get a solution vector
     * @param v
     */
    virtual void getSolVec(GlobalVectorType &v) = 0;

    /**
     * set a solution vector
     * @param v
     */
    virtual void setSolVec(const GlobalVectorType &v) = 0;

    /**
     * add RHS
     * @param dofId
     * @param list_discrete_pt_id
     * @param list_rhs_values
     * @param fac
     */
    virtual void addRHSVec(
                    std::size_t dofId,
                    const std::vector<std::size_t> &list_discrete_pt_id,
                    const std::vector<double> &list_rhs_values,
                    double fac) = 0;

    /**
     * add RHS vector
     * @param v     RHS vector
     * @param fac   Signed scaling factor
     */
    virtual void addRHSVec(const GlobalVectorType &v, double fac) = 0;
};


} //DiscreteLib

#endif //IDISCRETELINEARSYSTEM_H_
