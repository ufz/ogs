/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-17
 * \brief  Definition of the GEOTYPE enumeration.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOTYPE_H_
#define GEOTYPE_H_

#include <string>

namespace GeoLib {

/**
 * \ingroup GeoLib
 */

enum class GEOTYPE {
    POINT,
    POLYLINE,
    SURFACE
};

std::string convertGeoTypeToString (GEOTYPE geo_type);

} // end namespace GeoLib

#endif /* GEOTYPE_H_ */
