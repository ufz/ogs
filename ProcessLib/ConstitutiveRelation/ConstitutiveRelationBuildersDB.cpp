/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstitutiveRelationBuildersDB.h"

#include "IdealGasBinaryMixtureDensity.h"

namespace ProcessLib {
namespace ConstitutiveRelation {

ConstitutiveRelationBuildersDB&
ConstitutiveRelationBuildersDB::singleton()
{
    // cf. http://stackoverflow.com/a/3746390
    static ConstitutiveRelationBuildersDB instance;
    return instance;
}

#define OGS_REGISTER_CONSTITUTIVE_RELATION(CTIVE_REL, RETURN_TYPE, ...) \
    add(#CTIVE_REL, \
        std::unique_ptr<ConstitutiveRelationBuilder< \
            RETURN_TYPE, __VA_ARGS__>>(new CTIVE_REL ## Builder))

void ConstitutiveRelationBuildersDB::initialize()
{
    DBUG("initializing constitutive relations database");

    OGS_REGISTER_CONSTITUTIVE_RELATION(IdealGasBinaryMixtureDensity,
                                       double, double, double, double);
}

} // namespace ConstitutiveRelation
} // namespace ProcessLib
