/**
 * \author Norbert Grunwald
 * \date   12.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef MATERIALLIB_MPL_COMPONENTS_CSALT_H_
#define MATERIALLIB_MPL_COMPONENTS_CSALT_H_

#include "../mpComponent.h"

namespace MaterialPropertyLib
{
/**
 * \class Salt
 * \brief A class for Salt derived from Component
 * \details This class can holds material constants and default
 * properties of ordinary salt (NaCl in this case).
 */
class Salt : public Component
{
public:
    Salt();
};

}  // MaterialPropertyLib

#endif /* MATERIALLIB_MPL_COMPONENTS_CSALT_H_ */
