/**
 * \file   SerialDiscreteLinearSystem.tpp
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
#ifndef SERIALDISCRETELINEARSYSTEM_TPP_
#define SERIALDISCRETELINEARSYSTEM_TPP_

#include <cassert>

#include "MathLib/LinAlg/Sparse/Sparsity.h"

#include "MeshLib/Mesh.h"

#include "DiscreteLib/Core/IDiscreteLinearSystem.h"
#include "DiscreteLib/Core/IDiscreteLinearSystemAssembler.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"


namespace DiscreteLib
{

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::SerialDiscreteLinearSystem(const MeshLib::Mesh* msh,  const DofEquationIdTable* dofMap, DiscreteResourceManager* dis_resource)
: AbstractMeshBasedDiscreteLinearSystem(msh, dofMap), _dis_resource(dis_resource)
{
    assert(dofMap->getNumberOfVariables()>0);
    MathLib::RowMajorSparsity* sparse = new MathLib::RowMajorSparsity();
    MySparsityBuilder sp_builder(getMesh(), *dofMap, *sparse);
    AbstractMeshBasedDiscreteLinearSystem::setSparsity(sparse);
    _linear_sys = new MyLinearSolverType(dofMap->getTotalNumberOfActiveDoFs(), sparse);
};

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::setZero()
{
    _linear_sys->setZero();
    _list_prescribed_eqs_id.clear();
    _list_prescribed_values.clear();
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::construct(IDiscreteLinearSystemAssembler& assemler)
{
    assert(getDofEquationIdTable().getNumberOfVariables()>0);

    assemler.assembly(getMesh(), AbstractMeshBasedDiscreteLinearSystem::getDofEquationIdTable(), *_linear_sys);

    //apply 1st bc
    _linear_sys->setKnownSolution(_list_prescribed_eqs_id, _list_prescribed_values);
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::setKnownSolution(std::size_t varId, const std::vector<std::size_t> &list_discrete_pt_id, const std::vector<double> &list_prescribed_values)
{
    const std::size_t n = list_discrete_pt_id.size();
    const DofEquationIdTable &dofmap = getDofEquationIdTable();
    const std::size_t msh_id = getMesh().getID();
    for (std::size_t i=0; i<n; i++) {
        const DoF dof(varId, msh_id, list_discrete_pt_id[i]);
        if (dofmap.isActiveDoF(dof)) {
            _list_prescribed_eqs_id.push_back(dofmap.mapEqsID(dof));
            _list_prescribed_values.push_back(list_prescribed_values[i]);
        }
    }
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::addRHSVec(std::size_t varId, const std::vector<std::size_t> &list_discrete_pt_id, const std::vector<double> &list_rhs_values, double fkt)
{
    const std::size_t n = list_discrete_pt_id.size();
    const DofEquationIdTable &dofmap = getDofEquationIdTable();
    const std::size_t msh_id = getMesh().getID();
    for (std::size_t i=0; i<n; i++) {
        const DoF dof(varId, msh_id, list_discrete_pt_id[i]);
        if (dofmap.isActiveDoF(dof)) {
            _linear_sys->addRHSVec(dofmap.mapEqsID(dof), list_rhs_values[i]*fkt);
        }
    }
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::addRHSVec(const GlobalVectorType &v, double fkt)
{
    const std::size_t n = v.size();
    for (std::size_t i=0; i<n; i++) {
        _linear_sys->addRHSVec(i, v[i]*fkt);
    }
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::solve()
{
    _linear_sys->solve();
}


template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::getSolVec(std::vector<double> &x)
{
    x.resize(_linear_sys->getDimension());
    for (std::size_t i=0; i<x.size(); i++)
        x[i] = _linear_sys->getSolVec(i);
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::getSolVec(GlobalVectorType &x)
{
    for (std::size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
        x[i] = _linear_sys->getSolVec(i);
};

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::setSolVec(const GlobalVectorType &x)
{
    for (std::size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
        _linear_sys->setSolVec(i, x[i]);
};

} //end

#endif //SERIALDISCRETELINEARSYSTEM_TPP_
