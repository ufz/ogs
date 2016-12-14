/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowMaterialProperties.h
 *
 * Created on August 18, 2016, 11:03 AM
 */

#pragma once

#include <memory>
#include <Eigen/Dense>

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"

namespace MaterialLib
{
namespace Fluid
{
class FluidProperties;
}
}

namespace MaterialLib
{
namespace PorousMedium
{
class Porosity;
class Storage;
}
}

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
class SpatialPosition;

namespace LiquidFlow
{
class LiquidFlowMaterialProperties
{
public:
    typedef MaterialLib::Fluid::FluidProperty::ArrayType ArrayType;

    LiquidFlowMaterialProperties(
        BaseLib::ConfigTree const& config,
        bool const has_material_ids,
        MeshLib::PropertyVector<int> const& material_ids);

    void setMaterialID(const SpatialPosition& pos);

    /**
     * \brief Compute the coefficient of the mass term by
     *      \f[
     *           n \frac{partial \rho_l}{\partial p} + \beta_s
     *      \f]
     *     where \f$n\f$ is the porosity, \f$rho_l\f$ is the liquid density,
     *     \f$bata_s\f$ is the storage.
     * \param t                  Time.
     * \param pos                Position of element.
     * \param p                  Pressure value.
     * \param T                  Temperature value.
     * \param porosity_variable  The first variable for porosity model, and it
     *                           passes a double type value that could be
     *                           saturation, and invariant of stress or strain.
     * \param storage_variable   Variable for storage model.
     */
    double getMassCoefficient(const double t, const SpatialPosition& pos,
                              const double p, const double T,
                              const double porosity_variable,
                              const double storage_variable) const;

    Eigen::MatrixXd const& getPermeability(const double t,
                                           const SpatialPosition& pos,
                                           const int dim) const;

    double getLiquidDensity(const double p, const double T) const;

    double getViscosity(const double p, const double T) const;

private:
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> _fluid_properties;

    /// A flag to indicate whether the reference member, _material_ids,
    /// is not assigned.
    const bool _has_material_ids;
    /** Use porous medium models for different material zones.
     *  Material IDs must be given as mesh element properties.
     */
    MeshLib::PropertyVector<int> const& _material_ids;

    int _current_material_id = 0;
    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;

    // Note: For the statistical data of porous media, they could be read from
    // vtu files directly. This can be done by using property vectors directly.
    // Such property vectors will be added here if they are needed.
};

}  // end of namespace
}  // end of namespace
