/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "IdealGasBinaryMixtureDensity.h"

#include "BaseLib/ConfigTree.h"

namespace ProcessLib {
namespace ConstitutiveRelation {

std::unique_ptr<ConstitutiveRelation<double, double, double, double>>
IdealGasBinaryMixtureDensityBuilder::
createConstitutiveRelation(BaseLib::ConfigTree const& config)
{
    (void) config;
    return nullptr;
}

} // namespace ConstitutiveRelation
} // namespace ProcessLib
