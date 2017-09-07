/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef MATERIALLIB_MPL_MPPHASE_H_
#define MATERIALLIB_MPL_MPPHASE_H_

#include "mpComponent.h"

#include <string>
#include <vector>

namespace MaterialPropertyLib
{
class Phase
{
private:
    std::string _name;
    std::vector<Component*> _components;
public:
    Phase();
};

} //MaterialPropertyLib



#endif /* MATERIALLIB_MPL_MPPHASE_H_ */
