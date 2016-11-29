/**
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CompositeDensityViscosityModel.cpp
 *
 * Created on November 29, 2016, 3:19 PM
 */

#include "CompositeDensityViscosityModel.h"

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
CompositeDensityViscosityModel::CompositeDensityViscosityModel(
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& viscosity)
    : _density(std::move(density)), _viscosity(std::move(viscosity))
{
}

inline double CompositeDensityViscosityModel::getValue(
    const PropertyType property_type, const ArrayType& variable_values) const
{
    switch (property_type)
    {
        case PropertyType::Density:
            return _density->getValue(variable_values);
        case PropertyType::Vicosity:
            return _viscosity->getValue(variable_values);
        default:
            return 0.;
    }
}

inline double CompositeDensityViscosityModel::getdValue(
    const PropertyType property_type,
    const ArrayType& variable_values,
    const PropertyVariableType variable_type) const
{
    switch (property_type)
    {
        case PropertyType::Density:
            return _density->getdValue(variable_values, variable_type);
        case PropertyType::Vicosity:
            return _viscosity->getdValue(variable_values, variable_type);
        default:
            return 0.;
    }
}

}  // end namespace
}  // end namespace
