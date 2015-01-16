/**
 * \file
 * \date   2015-01-16
 * \brief  Definition of the MathPoint class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHPOINT_H_
#define MATHPOINT_H_

#include "TemplatePoint.h"

namespace MathLib
{
typedef MathLib::TemplatePoint<double,3> MathPoint;
} // end namespace MathLib

#endif /* MATHPOINT_H_ */

