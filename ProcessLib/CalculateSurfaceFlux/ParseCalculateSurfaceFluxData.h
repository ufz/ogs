
/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
void parseCalculateSurfaceFluxData(BaseLib::ConfigTree const& config,
                                   std::string& mesh_name,
                                   std::string& property_name,
                                   std::string& out_fname);
}  // namespace ProcessLib
