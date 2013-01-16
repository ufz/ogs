
/**
 * \file   IElemenetWiseLinearSystemLocalAssembler.h
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

#ifndef IELEMENETWISELINEARSYSTEMLOCALASSEMBLER_H_
#define IELEMENETWISELINEARSYSTEMLOCALASSEMBLER_H_

#include "DiscreteLib/DiscreteEnums.h"

namespace MeshLib
{
class Element;
}

namespace DiscreteLib
{

/**
 * \brief Interface of element-wise local assembler classes
 */
class IElemenetWiseLinearSystemLocalAssembler
{
public:
    ///
    virtual ~IElemenetWiseLinearSystemLocalAssembler(){};

    /**
     * assemble a local linear equation for the given element
     * @param e     Mesh element
     * @param eqs   Local linear system, i.e. Ax=b
     */
    virtual void assembly(const MeshLib::Element &e, LocalLinearSystem &eqs) = 0;
};

} //end

#endif //IELEMENETWISELINEARSYSTEMLOCALASSEMBLER_H_
