/**
 * \file   SequentialElementWiseLinearSystemAssembler.h
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
#ifndef SEQUENTIALELEMENTWISELINEARSYSTEMASSEMBLER_H_
#define SEQUENTIALELEMENTWISELINEARSYSTEMASSEMBLER_H_

#include "MeshLib/Mesh.h"
#include "DiscreteLib/Core/IDiscreteLinearSystemAssembler.h"
#include "DiscreteLib/ElementWiseManipulator/IElemenetWiseLinearSystemLocalAssembler.h"
#include "DiscreteLib/DoF/DofEquationIdTable.h"

namespace DiscreteLib
{

/**
 * \brief Element-based linear system assembler
 *
 * \tparam T_UPDATER    Local updater class
 * \tparam T_SOLVER     Linear system class
 */
template <class T_LOCAL_ASSEMBLER, class T_SOLVER>
class SequentialElementWiseLinearSystemAssembler : public IDiscreteLinearSystemAssembler
{
public:
    typedef T_LOCAL_ASSEMBLER LocalAssembler;
    typedef T_SOLVER SolverType;

    /**
     * Constructor
     * @param updater
     */
    explicit SequentialElementWiseLinearSystemAssembler(const LocalAssembler &updater) : _e_assembler(updater) {};

    ///
    virtual ~SequentialElementWiseLinearSystemAssembler(){};

    /**
     * Conduct the element by element assembly procedure
     * @param msh
     * @param dofEquationIdTable
     * @param eqs
     */
    virtual void assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, MathLib::ISystemOfLinearEquations &eqs)
    {
        //TODO this function is not needed but how interface?
        assembly(msh, dofEquationIdTable, *((SolverType*)&eqs));
    }

    /**
     * Conduct the element by element assembly procedure
     * @param msh                   Mesh
     * @param dofEquationIdTable    DoF mapping
     * @param eqs                   Linear system
     */
    void assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, SolverType &eqs);

private:
    LocalAssembler _e_assembler;
};

template <class T1, class T2>
void SequentialElementWiseLinearSystemAssembler<T1, T2>::assembly(const MeshLib::Mesh &msh, const DofEquationIdTable &dofEquationIdTable, SolverType &eqs)
{
    const std::size_t n_ele = msh.getNElements();
    for (std::size_t i=0; i<n_ele; i++) {
        const MeshLib::Element *e = msh.getElement(i);
        _e_assembler.update(*e, dofEquationIdTable, eqs);
    }
};

}

#endif //SEQUENTIALELEMENTWISELINEARSYSTEMASSEMBLER_H_
