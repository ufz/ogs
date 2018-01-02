/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowMaterialProperties.cpp
 *
 * Created on August 18, 2016, 11:49 AM
 */

#include "LiquidFlowMaterialProperties.h"

#include <logog/include/logog.hpp>

#include "BaseLib/reorderVector.h"

#include "MeshLib/PropertyVector.h"

#include "ProcessLib/Parameter/SpatialPosition.h"

#include "MaterialLib/PorousMedium/Permeability/Permeability.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"

namespace ProcessLib
{
namespace LiquidFlow
{
int LiquidFlowMaterialProperties::getMaterialID(
    const SpatialPosition& pos) const
{
    if (!_has_material_ids)
    {
        return 0;
    }

    assert(pos.getElementID().get() < _material_ids.size());
    return _material_ids[pos.getElementID().get()];
}

double LiquidFlowMaterialProperties::getLiquidDensity(const double p,
                                                      const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _fluid_properties->getValue(
        MaterialLib::Fluid::FluidPropertyType::Density, vars);
}

double LiquidFlowMaterialProperties::getViscosity(const double p,
                                                  const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _fluid_properties->getValue(
        MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
}

double LiquidFlowMaterialProperties::getPorosity(const int material_id,
                                                 const double t,
                                                 const SpatialPosition& pos,
                                                 const double porosity_variable,
                                                 const double T) const
{
    return _porosity_models[material_id]->getValue(t, pos, porosity_variable,
                                                   T);
}

double LiquidFlowMaterialProperties::getMassCoefficient(
    const int material_id, const double t, const SpatialPosition& pos,
    const double p, const double T, const double porosity_variable,
    const double storage_variable) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    const double drho_dp = _fluid_properties->getdValue(
        MaterialLib::Fluid::FluidPropertyType::Density, vars,
        MaterialLib::Fluid::PropertyVariableType::p);
    const double rho = _fluid_properties->getValue(
        MaterialLib::Fluid::FluidPropertyType::Density, vars);
    assert(rho > 0.);

    const double porosity =
        _porosity_models[material_id]->getValue(t, pos, porosity_variable, T);
    const double storage =
        _storage_models[material_id]->getValue(storage_variable);
    return porosity * drho_dp / rho + storage;
}

Eigen::MatrixXd const& LiquidFlowMaterialProperties::getPermeability(
    const int material_id, const double t, const SpatialPosition& pos,
    const int /*dim*/) const
{
    return _intrinsic_permeability_models[material_id]->getValue(t, pos, 0.0,
                                                                 0.0);
}

}  // end of namespace
}  // end of namespace
