// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

#include "MeshLib/PropertyVector.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct DensityBasedSurfaceSite
{
    DensityBasedSurfaceSite(std::string name_,
                            double const site_density_,
                            double const specific_surface_area_,
                            double const mass_)
        : name(std::move(name_)),
          site_density(site_density_),
          specific_surface_area(specific_surface_area_),
          mass(mass_)
    {
    }

    std::string const name;
    double const site_density;
    double const specific_surface_area;
    double const mass;
};

struct MoleBasedSurfaceSite
{
    MoleBasedSurfaceSite(std::string name_,
                         MeshLib::PropertyVector<double>* const molality_)
        : name(std::move(name_)), molality(molality_)
    {
    }

    std::string const name;
    MeshLib::PropertyVector<double>* molality;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
