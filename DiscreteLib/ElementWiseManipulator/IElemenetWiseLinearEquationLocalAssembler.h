
/**
 * \file   IElemenetWiseLinearSystemLocalAssembler
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

#ifndef IELEMENETWISELINEAREQUATIONLOCALASSEMBLER_H_
#define IELEMENETWISELINEAREQUATIONLOCALASSEMBLER_H_

#include "DiscreteLib/DiscreteEnums.h"

namespace MeshLib
{
class Element;
}

namespace DiscreteLib
{

/**
 * \brief Interface of all element local assembler classes
 */
class IElemenetWiseLinearSystemLocalAssembler
{
public:
    virtual ~IElemenetWiseLinearSystemLocalAssembler(){};

    /// assemble a local linear equation for the given element
    virtual void assembly(const MeshLib::Element &e, LocalLinearSystem &eqs) = 0;
};

} //end

#endif //IELEMENETWISELINEAREQUATIONLOCALASSEMBLER_H_
