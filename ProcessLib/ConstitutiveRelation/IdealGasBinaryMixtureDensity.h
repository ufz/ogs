/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CONSTITUTIVERELATION_IDEALGASBINARYMIXTUREDENSITY_H
#define PROCESSLIB_CONSTITUTIVERELATION_IDEALGASBINARYMIXTUREDENSITY_H

#include "ConstitutiveRelationBuilder.h"

namespace ProcessLib
{
namespace ConstitutiveRelation
{

OGS_DECLARE_CONSTITUTIVE_RELATION_BUILDER(IdealGasBinaryMixtureDensity, double,
                                          double, double, double);

}  // namespace ConstitutiveRelation
}  // namespace ProcessLib

#endif  // PROCESSLIB_CONSTITUTIVERELATION_IDEALGASBINARYMIXTUREDENSITY_H
