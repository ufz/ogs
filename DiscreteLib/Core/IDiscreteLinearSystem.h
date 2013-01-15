/**
 * \file   IDiscreteLinearSystem.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
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
 * \brief Interface for discrete linear equations
 */
class IDiscreteLinearSystem: public DiscreteObjectWithID
{
public:
    typedef IDiscreteVector<double> GlobalVectorType;

    ///
    virtual ~IDiscreteLinearSystem(){};

    /**
     * initialize a linear system
     */
    virtual void initialize() = 0;

    /**
     * construct this linear system
     *
     * @param assembler
     */
    virtual void construct(IDiscreteLinearSystemAssembler& assembler) = 0;

    /**
     *
     */
    virtual void solve() = 0;

    /**
     *
     * @param x
     */
    virtual void getGlobalX(std::vector<double> &x) = 0;

    /**
     *
     * @return
     */
    virtual double* getLocalX() = 0;

    /**
     *
     * @param v
     */
    virtual void getX(GlobalVectorType &v) = 0;

    /**
     *
     * @param v
     */
    virtual void setX(const GlobalVectorType &v) = 0;

    /**
     *
     * @return
     */
    virtual DofEquationIdTable* getDofMapManger() const = 0;

    /**
     *
     * @param dofId
     * @param list_discrete_pt_id
     * @param list_prescribed_values
     */
    virtual void setPrescribedDoF(std::size_t dofId, std::vector<std::size_t> &list_discrete_pt_id, std::vector<double> &list_prescribed_values) = 0;

    /**
     *
     * @param dofId
     * @param list_discrete_pt_id
     * @param list_rhs_values
     * @param fkt
     */
    virtual void addRHS(std::size_t dofId, std::vector<std::size_t> &list_discrete_pt_id, std::vector<double> &list_rhs_values, double fkt) = 0;

    /**
     *
     * @param v
     * @param fkt
     */
    virtual void addRHS(const GlobalVectorType &v, double fkt) = 0;
};


} //DiscreteLib

#endif //IDISCRETELINEAREQUATION_H_
