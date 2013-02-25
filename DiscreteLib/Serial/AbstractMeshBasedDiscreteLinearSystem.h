/**
 * \file   AbstractMeshBasedDiscreteLinearSystem.h
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
#ifndef ABSTRACTMESHBASEDDISCRETELINEARSYSTEM_H_
#define ABSTRACTMESHBASEDDISCRETELINEARSYSTEM_H_

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/Sparse/Sparsity.h"
#include "DiscreteLib/Core/IDiscreteLinearSystem.h"

namespace MeshLib
{
class Mesh;
}

namespace DiscreteLib
{
class DofEquationIdTable;

/**
 * \brief Abstract class for a single mesh-based discrete linear equations
 * 
 * - Mesh
 * - DoF manager
 * - Sparse pattern
 */
class AbstractMeshBasedDiscreteLinearSystem : public IDiscreteLinearSystem
{
public:
    /**
     * Constructor
     * @param msh
     * @param dofManager
     */
    AbstractMeshBasedDiscreteLinearSystem(const MeshLib::Mesh* msh, const DofEquationIdTable* dofManager)
    : _msh(msh), _dofManager(dofManager), _sparsity(nullptr)
    {
    }

    /**
     *
     */
    virtual ~AbstractMeshBasedDiscreteLinearSystem()
    {
        if (_sparsity!=nullptr) {
            delete  _sparsity;
            _sparsity = nullptr;
        }
    }

    /// return a mesh
    const MeshLib::Mesh& getMesh() const { return *_msh; }

    /// return DoF map table
    virtual const DofEquationIdTable& getDofEquationIdTable() const { return *_dofManager; }

    /// return a sparse pattern of equation matrix
    const MathLib::RowMajorSparsity& getSparsity() const { return *_sparsity; };

protected:
    /// set sparse pattern
    void setSparsity(MathLib::RowMajorSparsity* sp) { _sparsity = sp; }

private:
    DISALLOW_COPY_AND_ASSIGN(AbstractMeshBasedDiscreteLinearSystem);

private:
    const MeshLib::Mesh* _msh;
    const DofEquationIdTable* _dofManager;
    MathLib::RowMajorSparsity* _sparsity;
};


} //end

#endif //ABSTRACTMESHBASEDDISCRETELINEARSYSTEM_H_
