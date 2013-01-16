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

#include "MathLib/LinAlg/Sparse/Sparsity.h"

#include "MeshLib/Mesh.h"

#include "DiscreteLib/Core/IDiscreteLinearSystem.h"
#include "DiscreteLib/Core/IDiscreteLinearSystemAssembler.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"


namespace DiscreteLib
{

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>* SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::createInstance(IDiscreteSystem &dis_sys, const MeshLib::Mesh* msh, DofEquationIdTable* dofManager)
{
    MySerialDiscreteLinearSystem *eqs = new MySerialDiscreteLinearSystem(msh, dofManager);
    dis_sys.addLinearSystem(eqs);
    return eqs;
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::SerialDiscreteLinearSystem(const MeshLib::Mesh* msh,  DofEquationIdTable* dofMap)
    : AbstractMeshBasedDiscreteLinearSystem(msh, dofMap)
{
    assert(dofMap->getNumberOfVariables()>0);
    MathLib::RowMajorSparsity* sparse = new MathLib::RowMajorSparsity();
    MySparsityBuilder sp_builder(*getMesh(), *dofMap, *sparse);
    AbstractMeshBasedDiscreteLinearSystem::setSparsity(sparse);
    _eqs = new MyLinearSolverType(dofMap->getTotalNumberOfActiveDoFs(), sparse);
};

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::setZero()
{
    _eqs->setZero();
    _list_prescribed_dof_id.clear();
    _list_prescribed_values.clear();
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::construct(IDiscreteLinearSystemAssembler& assemler)
{
    assert(getDofEquationIdTable()->getNumberOfVariables()>0);

    assemler.assembly(*getMesh(), *AbstractMeshBasedDiscreteLinearSystem::getDofEquationIdTable(), *_eqs);

    //apply 1st bc
    _eqs->setKnownSolution(_list_prescribed_dof_id, _list_prescribed_values);
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::setKnownSolution(size_t varId, const std::vector<size_t> &list_discrete_pt_id, const std::vector<double> &list_prescribed_values)
{
    const size_t n = list_discrete_pt_id.size();
    const DofEquationIdTable* dofmap = getDofEquationIdTable();
    const size_t msh_id = getMesh()->getID();
    for (size_t i=0; i<n; i++) {
        size_t pt_id = list_discrete_pt_id[i];
        if (dofmap->isActiveDoF(varId, msh_id, pt_id)) {
            _list_prescribed_dof_id.push_back(dofmap->mapEqsID(varId, msh_id, pt_id));
            _list_prescribed_values.push_back(list_prescribed_values[i]);
        }
    }
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::addRHSVec(size_t dofId, const std::vector<size_t> &list_discrete_pt_id, const std::vector<double> &list_rhs_values, double fkt)
{
    const size_t n = list_discrete_pt_id.size();
    const DofEquationIdTable* dofmap = getDofEquationIdTable();
    const size_t msh_id = getMesh()->getID();
    for (size_t i=0; i<n; i++) {
        size_t pt_id = list_discrete_pt_id[i];
        if (dofmap->isActiveDoF(dofId, msh_id, pt_id)) {
            _eqs->addRHSVec(dofmap->mapEqsID(dofId, msh_id, pt_id), list_rhs_values[i]*fkt);
        }
    }
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::addRHSVec(const GlobalVectorType &v, double fkt)
{
    const size_t n = v.size();
    for (size_t i=0; i<n; i++) {
        _eqs->addRHSVec(i, v[i]*fkt);
    }
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::solve()
{
    _eqs->solve();
}


template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::getSolVec(std::vector<double> &x)
{
    x.resize(_eqs->getDimension());
    for (size_t i=0; i<x.size(); i++)
        x[i] = _eqs->getSolVec(i);
}

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::getSolVec(GlobalVectorType &x)
{
    for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
        x[i] = _eqs->getSolVec(i);
};

template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
void SerialDiscreteLinearSystem<T_LINEAR_SOLVER,T_SPARSITY_BUILDER>::setSolVec(const GlobalVectorType &x)
{
    for (size_t i=x.getRangeBegin(); i<x.getRangeEnd(); i++)
        _eqs->setSolVec(i, x[i]);
};

} //end

#endif //SERIALDISCRETELINEARSYSTEM_TPP_
