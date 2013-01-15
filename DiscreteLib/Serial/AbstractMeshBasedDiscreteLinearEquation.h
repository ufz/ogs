/**
 * \file   AbstractMeshBasedDiscreteLinearEquation.h
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
#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Sparse/Sparsity.h"
#include "MeshLib/Mesh.h"
#include "DiscreteLib/Core/IDiscreteLinearSystem.h"

namespace DiscreteLib
{

/**
 * \brief Abstract class for mesh-based discrete linear equations
 * 
 * - Mesh
 * - DoF manager
 * - Sparse pattern
 */
class AbstractMeshBasedDiscreteLinearEquation : public IDiscreteLinearSystem
{
public:
    /// \param  msh
    /// \param  dofManager
    AbstractMeshBasedDiscreteLinearEquation(const MeshLib::Mesh* msh, DofEquationIdTable* dofManager)
    : _msh(msh), _dofManager(dofManager), _sparsity(nullptr)
    {
    }

    virtual ~AbstractMeshBasedDiscreteLinearEquation()
    {
        if (_sparsity!=nullptr) {
            delete  _sparsity;
            _sparsity = nullptr;
        }
    }

    /// return a mesh
    const MeshLib::Mesh* getMesh() const { return _msh; }

    /// return DoF map table
    virtual DofEquationIdTable* getDofMapManger() const { return _dofManager; }

    /// return a sparse pattern of equation matrix
    MathLib::RowMajorSparsity* getSparsity() const { return _sparsity; };

protected:
    void setSparsity(MathLib::RowMajorSparsity* sp) { _sparsity = sp; }

private:
    DISALLOW_COPY_AND_ASSIGN(AbstractMeshBasedDiscreteLinearEquation);

private:
    const MeshLib::Mesh* _msh;
    DofEquationIdTable* _dofManager;
    MathLib::RowMajorSparsity* _sparsity;
};


} //end
