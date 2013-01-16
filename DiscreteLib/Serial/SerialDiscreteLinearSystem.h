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
 * \brief Non-parallel version of a mesh based linear system
 *        combined with a specific linear solver
 *
 * \tparam T_LINEAR_SOLVER      Linear solver class
 * \tparam T_SPARSITY_BUILDER   Sparsity builder
 * 
 */
template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
class SerialDiscreteLinearSystem
 : public AbstractMeshBasedDiscreteLinearSystem
{
public:
    friend class SerialDiscreteSystem;

    typedef T_LINEAR_SOLVER MyLinearSolverType;
    typedef T_SPARSITY_BUILDER MySparsityBuilder;
    typedef SerialDiscreteLinearSystem<MyLinearSolverType,MySparsityBuilder>
                MySerialDiscreteLinearSystem;

    /// zero the system
    virtual void setZero();

    /// construct the linear equation
    virtual void construct(IDiscreteLinearSystemAssembler& assemler);

    /// set prescribed dof
    virtual void setKnownSolution(  std::size_t varId,
                                    const std::vector<std::size_t> &list_discrete_pt_id,
                                    const std::vector<double> &list_prescribed_values);
   
    /// set additional RHS values
    virtual void addRHSVec( std::size_t dofId,
                            const std::vector<std::size_t> &list_discrete_pt_id,
                            const std::vector<double> &list_rhs_values,
                            double fkt=1.0);

    /// add RHS
    virtual void addRHSVec(const GlobalVectorType &v, double fkt=1.0);

    /// solve
    virtual void solve();

    /// get a solution vector
    virtual void getSolVec(std::vector<double> &x);

    /// get a solution vector
    virtual void getSolVec(GlobalVectorType &x);

    /// set a solution vector
    virtual void setSolVec(const GlobalVectorType &x);

    /// return a linear solver object
    MyLinearSolverType* getLinearSolver() {return _linear_sys;};

private:
    /**
     * constructor
     * @param msh           Mesh
     * @param dofMap        DoF mapping table
     * @param dis_resource  Discrete data resource manager
     */
    SerialDiscreteLinearSystem( const MeshLib::Mesh* msh,
                                const DofEquationIdTable* dofMap,
                                DiscreteResourceManager* dis_resource);

private:
    DISALLOW_COPY_AND_ASSIGN(SerialDiscreteLinearSystem);

private:
    DiscreteResourceManager* _dis_resource;
    MyLinearSolverType* _linear_sys;
    std::vector<std::size_t> _list_prescribed_eqs_id;
    std::vector<double> _list_prescribed_values;

};

} //end

#include "SerialDiscreteLinearSystem.tpp"

#endif //SERIALDISCRETELINEARSYSTEM_H_
