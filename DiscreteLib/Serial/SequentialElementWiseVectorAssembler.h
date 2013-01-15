/**
 * \file   SequentialElementWiseVectorAssembler.h
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

#include <vector>

#include "MeshLib/Mesh.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"
#include "DiscreteLib/Core/IDiscreteVectorAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseVectorLocalAssembler.h"

namespace DiscreteLib
{

/**
 * \brief Element-based discrete vector assembler classes
 */
template <class T_VALUE, class T_UPDATER>
class SequentialElementWiseVectorAssembler : public IDiscreteVectorAssembler<T_VALUE>
{
public:
    typedef typename IDiscreteVectorAssembler<T_VALUE>::VectorType GlobalVectorType;
    typedef T_UPDATER UpdaterType;

    ///
    SequentialElementWiseVectorAssembler(UpdaterType* a) : _e_assembler(a) {};
    virtual ~SequentialElementWiseVectorAssembler() {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param vec Discrete vector
    virtual void assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, GlobalVectorType &globalVec);

private:
    UpdaterType* _e_assembler;
};


template<class T1, class T2>
void SequentialElementWiseVectorAssembler<T1,T2>::assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, GlobalVectorType &globalVec)
{
    const std::size_t n_ele = msh.getNElements();

    for (std::size_t i=0; i<n_ele; i++) {
        MeshLib::Element *e = msh.getElement(i);
        _e_assembler->update(*e, dofEquationIdTable, globalVec);
    }
};

}
