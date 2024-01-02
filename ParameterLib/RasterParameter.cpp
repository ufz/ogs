/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RasterParameter.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

namespace ParameterLib
{
std::unique_ptr<ParameterBase> createRasterParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    std::vector<GeoLib::NamedRaster> const& named_rasters)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "Raster");

    auto const& named_raster = BaseLib::getIfOrError(
        named_rasters,
        [&name](auto const& named_raster)
        { return name == named_raster.raster_name; },
        "Could not find raster '" + name);

    DBUG("Using the raster '{}' for the raster parameter.", name);

    return std::make_unique<RasterParameter>(name, named_raster);
}

}  // namespace ParameterLib
