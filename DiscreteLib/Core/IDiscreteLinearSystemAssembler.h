/**
 * \file   IDiscreteLinearSystemAssembler.h
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

#ifndef IDISCRETELINEARSYSTEMASSEMBLER_H_
#define IDISCRETELINEARSYSTEMASSEMBLER_H_

namespace MeshLib 
{
class Mesh;
}

namespace MathLib
{
class ISystemOfLinearEquations;
}

namespace DiscreteLib
{

class DofEquationIdTable;

/**
 * \brief Interface of discrete system assembler classes
 */
class IDiscreteLinearSystemAssembler
{
public:
    ///
    virtual ~IDiscreteLinearSystemAssembler() {};

    /**
     * assemble a linear system
     *
     * @param msh                   Mesh object
     * @param dofEquationIdTable    DoF mapping table
     * @param eqs                   Linear system to be updated
     */
    virtual void assembly(
            const MeshLib::Mesh &msh,
            const DofEquationIdTable &dofEquationIdTable,
            MathLib::ISystemOfLinearEquations &eqs) = 0;
};

}

#endif //IDISCRETELINEARSYSTEMASSEMBLER_H_

