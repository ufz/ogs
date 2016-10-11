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

#ifndef OGS_LIQUIDFLOWMATERIALPROPERTIES_H
#define OGS_LIQUIDFLOWMATERIALPROPERTIES_H

#include <memory>

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"

#include "MeshLib/PropertyVector.h"

namespace MeshLib
{
class Mesh;
}

namespace MaterialLib
{
namespace PorousMedium
{
class Porosity;
class Storage;
}
}

namespace ProcessLib
{
namespace LiquidFlow
{
struct LiquidFlowMaterialProperties
{
    typedef MaterialLib::Fluid::FluidProperty::ArrayType ArrayType;

    explicit LiquidFlowMaterialProperties(BaseLib::ConfigTree const& config);

    /**
     * \brief Compute the coefficient of the mass term by
     *      \f[
     *           n \frac{partial \rho_l}{\partial p} + \beta_s
     *      \f]
     *     where \f$n\f$ is the porosity, \f$rho_l\f$ is the liquid density,
     *     \f$bata_s\f$ is the storage.
     * \param porosity_variable  The first variable for porosity model, and it
     *                           passes a double type value that could be
     *                           saturation, and invariant of stress or strain.
     * \param storage_variable   Variable for storage model.
     * \param p                  Pressure value
     * \param T                  Temperature value
     * \param material_group_id  Material ID of the element
     */
    double getMassCoefficient(const double p, const double T,
                              const double porosity_variable,
                              const double storage_variable,
                              const unsigned material_group_id = 0) const;

    double getLiquidDensity(const double p, const double T) const;

    double getViscosity(const double p, const double T) const;

    std::unique_ptr<MaterialLib::Fluid::FluidProperty> liquid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> viscosity;

    /// Porous medium properties of different material zones.
    /// The vector is left empty if the property data are given in vtu file,
    /// e.g for heterogeneous medium.
    std::vector<Eigen::MatrixXd> intrinsic_permeability;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>> porosity;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>> storage;

    // TODO: heterogeneous medium.
};

}  // end of namespace
}  // end of namespace
#endif /* LIQUIDFLOWMATERIALPROPERTIES_H */
