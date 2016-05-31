/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstitutiveRelationsDB.h"

namespace ProcessLib {
namespace ConstitutiveRelation {

bool db_dummy = false;

ConstitutiveRelationsDB&
ConstitutiveRelationsDB::singleton()
{
    // cf. http://stackoverflow.com/a/3746390
    static ConstitutiveRelationsDB instance;
    return instance;
}

} // namespace ConstitutiveRelation
} // namespace ProcessLib
