/**
 * \file   SequentialElementWiseLinearEquationAssembler.h
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

#include "MeshLib/Mesh.h"
#include "DiscreteLib/Core/IDiscreteLinearSystemAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseLinearEquationLocalAssembler.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"

namespace DiscreteLib
{

/**
 * \brief Element-based linear equation assembler
 */
template <class T_UPDATER, class T_SOLVER>
class SequentialElementWiseLinearEquationAssembler : public IDiscreteLinearSystemAssembler
{
public:
    typedef T_UPDATER UpdaterType;
    typedef T_SOLVER SolverType;

    /**
     *
     * @param updater
     */
    explicit SequentialElementWiseLinearEquationAssembler(UpdaterType* updater) : _e_assembler(updater) {};

    ///
    virtual ~SequentialElementWiseLinearEquationAssembler(){};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, SolverType &eqs);

    virtual void assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, MathLib::ISystemOfLinearEquations &eqs)
    {
        assembly(msh, dofEquationIdTable, *((SolverType*)&eqs));
    }

private:
    UpdaterType* _e_assembler;
};

template <class T1, class T2>
void SequentialElementWiseLinearEquationAssembler<T1, T2>::assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, SolverType &eqs)
{
    const std::size_t n_ele = msh.getNElements();
    for (std::size_t i=0; i<n_ele; i++) {
        MeshLib::Element *e = msh.getElement(i);
        _e_assembler->update(*e, dofEquationIdTable, eqs);
    }
};

}
