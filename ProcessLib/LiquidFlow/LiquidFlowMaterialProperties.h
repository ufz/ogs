/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowMaterialProperties.h
 *
 * Created on August 18, 2016, 11:03 AM
 */

#ifndef LIQUIDFLOWMATERIALPROPERTIES_H
#define LIQUIDFLOWMATERIALPROPERTIES_H

#include <memory>

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"

namespace ProcessLib
{
namespace LiquidFlow
{
struct LiquidFlowMaterialProperties
{
    LiquidFlowMaterialProperties(BaseLib::ConfigTree const& config);

    /**
     * \brief Compute the coefficient of the mass term by
     *      \f[
     *           n \frac{partial \rho_l}{\partial p} + \beta_s
     *      \f]
     *     where \f$n\f$ is the porosity, \f$rho_l\f$ is the liquid density,
     *     \f$bata_s\f$ is the storage.
     * @param p                  Pressure value
     * @param T                  Temperature value
     * @param material_group_id  Material ID of the element
     * @return
     */
    double getMassCoeffcient(const double p, const double T,
                             const unsigned material_group_id = 0);

    double getLiquidDensity(const double p, const double T)
    {
        vars[0] = T;
        vars[1] = p;
        return density_l->getValue(vars);
    }

    double getViscosity(const double p, const double T)
    {
        vars[0] = T;
        vars[1] = p;
        return viscosity->getValue(vars);
    }

    std::unique_ptr<MaterialLib::Fluid::FluidProperty> density_l;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> viscosity;

    /// Porous medium properties of different material zones.
    /// The vector is left empty if the property data are given in vtu file,
    /// e.g for heterogeneous medium.
    std::vector<MaterialLib::PorousMedium::CoefMatrix> intrinsic_permeabiliy;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>> porosity;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>> storage;

    double vars[2];
};

}  // end of namespace
}  // end of namespace
#endif /* LIQUIDFLOWMATERIALPROPERTIES_H */
