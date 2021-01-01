/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/optional/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "Surface.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<SurfaceSite> createSurface(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return {};
    }

    std::vector<SurfaceSite> surface;
    for (auto const& site_config :
         //! \ogs_file_param{prj__chemical_system__surface__site}
         config->getConfigSubtreeList("site"))
    {
        //! \ogs_file_param{prj__chemical_system__surface__site__name}
        auto name = site_config.getConfigParameter<std::string>("name");

        auto const site_density =
            //! \ogs_file_param{prj__chemical_system__surface__site__site_density}
            site_config.getConfigParameter<double>("site_density");

        auto const specific_surface_area =
            //! \ogs_file_param{prj__chemical_system__surface__site__specific_surface_area}
            site_config.getConfigParameter<double>("specific_surface_area");

        auto const mass =
            //! \ogs_file_param{prj__chemical_system__surface__site__mass}
            site_config.getConfigParameter<double>("mass");

        surface.emplace_back(
            std::move(name), site_density, specific_surface_area, mass);
    }

    return surface;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
