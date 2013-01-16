/**
 * \file   SerialDiscreteLinearSystem.h
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
#ifndef SERIALDISCRETELINEARSYSTEM_H_
#define SERIALDISCRETELINEARSYSTEM_H_

#include <vector>

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/SparsityBuilder/SparsityBuilderDummy.h"
#include "AbstractMeshBasedDiscreteLinearSystem.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{
class DofEquationIdTable;
class IDiscreteSystem;
class IDiscreteLinearSystemAssembler;

/**
 * \brief Implementation of mesh based linear equation classes 
 *        combined with a specific linear solver class
 *
 * \tparam T_LINEAR_SOLVER      Linear solver class
 * \tparam T_SPARSITY_BUILDER   Sparsity builder
 * 
 * - linear equation
 */
template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER=SparsityBuilderDummy>
class SerialDiscreteLinearSystem : public AbstractMeshBasedDiscreteLinearSystem
{
public:
    typedef T_LINEAR_SOLVER MyLinearSolverType;
    typedef T_SPARSITY_BUILDER MySparsityBuilder;
    typedef SerialDiscreteLinearSystem<MyLinearSolverType,MySparsityBuilder> MySerialDiscreteLinearSystem;

    /**
     * Create a SerialDiscreteLinearSystem object
     *
     * @param dis_sys           Discrete system object
     * @param msh               Mesh object
     * @param dofMap            DoF mapping table
     * @return
     */
    static MySerialDiscreteLinearSystem* createInstance(
                IDiscreteSystem &dis_sys,
                const MeshLib::Mesh* msh,
                DofEquationIdTable* dofMap);

    /// zero the system
    virtual void setZero();

    /// construct the linear equation
    virtual void construct(IDiscreteLinearSystemAssembler& assemler);

    /// set prescribed dof
    virtual void setKnownSolution(size_t varId, const std::vector<size_t> &list_discrete_pt_id, const std::vector<double> &list_prescribed_values);
   
    /// set additional RHS values
    virtual void addRHSVec(size_t dofId, const std::vector<size_t> &list_discrete_pt_id, const std::vector<double> &list_rhs_values, double fkt=1.0);

    /// add RHS
    virtual void addRHSVec(const GlobalVectorType &v, double fkt=1.0);

    /// solve
    virtual void solve();

    /// copy global vector of x
    virtual void getSolVec(std::vector<double> &x);

    /// get a global vector of x
    virtual void getSolVec(GlobalVectorType &x);

    ///
    virtual void setSolVec(const GlobalVectorType &x);

    /// return a linear solver object
    MyLinearSolverType* getLinearSolver() {return _eqs;};

private:
    /// constructor
    /// \param msh
    /// \param linear_solver
    /// \param dofManager
    SerialDiscreteLinearSystem(const MeshLib::Mesh* msh, DofEquationIdTable* dofMap);

private:
    DISALLOW_COPY_AND_ASSIGN(SerialDiscreteLinearSystem);

private:
    MyLinearSolverType* _eqs;
    std::vector<size_t> _list_prescribed_dof_id;
    std::vector<double> _list_prescribed_values;

};

} //end

#include "SerialDiscreteLinearSystem.tpp"

#endif //SERIALDISCRETELINEARSYSTEM_H_
