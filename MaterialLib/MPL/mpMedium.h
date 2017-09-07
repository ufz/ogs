/**
 * \author Norbert Grunwald
 * \date   07.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_MPL_MPMEDIUM_H_
#define MATERIALLIB_MPL_MPMEDIUM_H_

#include "BaseLib/ConfigTree.h"

namespace MaterialPropertyLib
{
class Medium
{
public:
    Medium(BaseLib::ConfigTree const&);
    void createPhases (BaseLib::ConfigTree const&);
    void createProperties (BaseLib::ConfigTree const&);
};

}

#endif /* MATERIALLIB_MPL_MPMEDIUM_H_ */
