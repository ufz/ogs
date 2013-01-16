/**
 * \file   IDiscreteVectorAssembler.h
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

#ifndef IDISCRETEVECTORASSEMBLER_H_
#define IDISCRETEVECTORASSEMBLER_H_

#include "IDiscreteVector.h"

namespace MeshLib 
{
class Mesh;
}

namespace DiscreteLib
{
class DofEquationIdTable;

/**
 * \brief Interface of discrete vector assembler classes
 *
 * \tparam T    Data type of a discrete vector, e.g. int, double
 */
template <class T>
class IDiscreteVectorAssembler
{
public:
    typedef IDiscreteVector<T> VectorType;

    virtual ~IDiscreteVectorAssembler() {};

    /**
     * Conduct the element by element assembly procedure
     *
     * @param msh                   Mesh object
     * @param dofEquationIdTable    DoF mapping table
     * @param vec                   Discrete vector to be updated
     */
    virtual void assembly(
            const MeshLib::Mesh &msh,
            const DofEquationIdTable &dofEquationIdTable,
            VectorType &vec) = 0;
};

}

#endif //IDISCRETEVECTORASSEMBLER_H_

