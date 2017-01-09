/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CREATEJACOBIANASSEMBLER_H
#define PROCESSLIB_CREATEJACOBIANASSEMBLER_H

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
class AbstractJacobianAssembler;

std::unique_ptr<AbstractJacobianAssembler> createJacobianAssembler(
    boost::optional<BaseLib::ConfigTree> const& config);
}  // ProcessLib

#endif  // PROCESSLIB_CREATEJACOBIANASSEMBLER_H
