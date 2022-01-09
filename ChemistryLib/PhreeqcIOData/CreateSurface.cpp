/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSurface.h"

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"
#include "Surface.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<std::variant<DensityBasedSurfaceSite, MoleBasedSurfaceSite>>
createSurface(std::optional<BaseLib::ConfigTree> const& config,
              MeshLib::Mesh& mesh)
{
    if (!config)
    {
        return {};
    }

    std::vector<std::variant<DensityBasedSurfaceSite, MoleBasedSurfaceSite>>
        surface;

    auto const surface_site_unit =
        //! \ogs_file_attr{prj__chemical_system__surface__site_unit}
        config->getConfigAttribute<std::string>("site_unit", "mole");

    if (surface_site_unit == "density")
    {
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

            surface.push_back(DensityBasedSurfaceSite(
                std::move(name), site_density, specific_surface_area, mass));
        }

        return surface;
    }

    if (surface_site_unit == "mole")
    {
        for (auto const& site_config :
             //! \ogs_file_param{prj__chemical_system__surface__site}
             config->getConfigSubtreeList("site"))
        {
            //! \ogs_file_param{prj__chemical_system__surface__site__name}
            auto name = site_config.getConfigParameter<std::string>("name");

            auto const molality = MeshLib::getOrCreateMeshProperty<double>(
                mesh, name, MeshLib::MeshItemType::IntegrationPoint, 1);

            surface.push_back(MoleBasedSurfaceSite(std::move(name), molality));
        }

        return surface;
    }

    OGS_FATAL("Surface site unit should be either of 'density' or 'mole'.");
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
