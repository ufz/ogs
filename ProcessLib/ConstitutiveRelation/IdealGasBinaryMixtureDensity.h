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

#include "ConstitutiveRelationBuildersDB.h"

namespace ProcessLib {
namespace ConstitutiveRelation {

// extern const AddConstitutiveRelationResult dummy;
// extern const bool dummy;

class IdealGasBinaryMixtureDensityBuilder final
        : public ConstitutiveRelationBuilder<double, double, double, double>
{
public:
    std::unique_ptr<ConstitutiveRelation<double, double, double, double>>
    createConstitutiveRelation(BaseLib::ConfigTree const& config) override;
};

} // namespace ConstitutiveRelation
} // namespace ProcessLib

#endif // PROCESSLIB_CONSTITUTIVERELATION_IDEALGASBINARYMIXTUREDENSITY_H
