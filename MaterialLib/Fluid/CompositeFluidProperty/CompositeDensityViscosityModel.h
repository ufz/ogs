/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CompositeDensityViscosityModel.h
 *
 * Created on November 29, 2016, 3:19 PM
 */

#ifndef OGS_COMPOSITE_DENSITY_VISCOSITY_MODEL_H
#define OGS_COMPOSITE_DENSITY_VISCOSITY_MODEL_H

#include <memory>

#include "MaterialLib/Fluid/CompositeFluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
class FluidProperty;

/// A class contains density and viscosity model.
class CompositeDensityViscosityModel final : public CompositeFluidProperty
{
public:
    CompositeDensityViscosityModel(
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity);

    /**
     *  Get the value of a Property.
     *  \param property_type   Property type.
     *  \param variable_values An array of variables. The order of its elements
     *                         is given in enum class PropertyVariableType.
     */
    double getValue(const PropertyType property_type,
                    const ArrayType& variable_values) const override;

    /**
     *  Get the partial differential of a property.
     *  \param property_type   Property type.
     *  \param variable_values An array of variables. The order of its elements
     *                         is given in enum class PropertyVariableType.
     *  \param variable_type   Variable type
     */
    double getdValue(const PropertyType property_type,
                     const ArrayType& variable_values,
                     const PropertyVariableType variable_type) const override;

private:
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _viscosity;
};

}  // end namespace
}  // end namespace
#endif /* OGS_COMPOSITE_DENSITY_VISCOSITY_MODEL_H */
