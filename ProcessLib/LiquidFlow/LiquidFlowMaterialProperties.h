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
#include <vector>

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"

#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

namespace MaterialLib
{
namespace Fluid
{
class FluidProperties;
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
template <typename T>
struct Parameter;

class SpatialPosition;

namespace LiquidFlow
{
/** This class has a collection of material properties of liquid flow process
 *  and it provides the computation of the coefficients of the PDE of liquid
 *  flow process.
 */
class LiquidFlowMaterialProperties
{
public:
    typedef MaterialLib::Fluid::FluidProperty::ArrayType ArrayType;

    LiquidFlowMaterialProperties(
        std::unique_ptr<MaterialLib::Fluid::FluidProperties>&& fluid_properties,
        std::vector<Eigen::MatrixXd>&& intrinsic_permeability_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>&&
            porosity_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>&&
            storage_models,
        bool const has_material_ids,
        MeshLib::PropertyVector<int> const& material_ids)
        : _has_material_ids(has_material_ids),
          _material_ids(material_ids),
          _fluid_properties(std::move(fluid_properties)),
          _intrinsic_permeability_models(
              std::move(intrinsic_permeability_models)),
          _porosity_models(std::move(porosity_models)),
          _storage_models(std::move(storage_models))
    {
    }

    int getMaterialID(const SpatialPosition& pos) const;

    /**
     * \brief Compute the coefficient of the mass term by
     *      \f[
     *           n \frac{partial \rho_l}{\partial p} + \beta_s
     *      \f]
     *     where \f$n\f$ is the porosity, \f$rho_l\f$ is the liquid density,
     *     \f$bata_s\f$ is the storage.
     * \param material_id        Material index.
     * \param t                  Time.
     * \param pos                Position of element.
     * \param p                  Pressure value.
     * \param T                  Temperature value.
     * \param porosity_variable  The first variable for porosity model, and it
     *                           passes a double type value that could be
     *                           saturation, and invariant of stress or strain.
     * \param storage_variable   Variable for storage model.
     */
    double getMassCoefficient(const int material_id, const double t,
                              const SpatialPosition& pos, const double p,
                              const double T, const double porosity_variable,
                              const double storage_variable) const;

    Eigen::MatrixXd const& getPermeability(const int material_id,
                                           const double t,
                                           const SpatialPosition& pos,
                                           const int dim) const;

    double getLiquidDensity(const double p, const double T) const;

    double getViscosity(const double p, const double T) const;

private:
    /// A flag to indicate whether the reference member, _material_ids,
    /// is not assigned.
    const bool _has_material_ids;
    /** Use porous medium models for different material zones.
     *  Material IDs must be given as mesh element properties.
     */
    MeshLib::PropertyVector<int> const& _material_ids;

    const std::unique_ptr<MaterialLib::Fluid::FluidProperties>
        _fluid_properties;

    const std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    const std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    const std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;

    // Note: For the statistical data of porous media, they could be read from
    // vtu files directly. This can be done by using property vectors directly.
    // Such property vectors will be added here if they are needed.
};

}  // end of namespace
}  // end of namespace
