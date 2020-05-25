/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TwoPhaseFlowWithPPMaterialProperties.h"

#include <boost/math/special_functions/pow.hpp>
#include <utility>
#include "BaseLib/Logging.h"

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowWithPPMaterialProperties::TwoPhaseFlowWithPPMaterialProperties(
    MeshLib::PropertyVector<int> const& material_ids,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& liquid_density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& liquid_viscosity,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& gas_density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>&& gas_viscosity,
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>&&
        intrinsic_permeability_models,
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>&&
        porosity_models,
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>&&
        storage_models,
    std::vector<std::unique_ptr<
        MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
        capillary_pressure_models,
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
        relative_permeability_models)
    : material_ids_(material_ids),
      liquid_density_(std::move(liquid_density)),
      liquid_viscosity_(std::move(liquid_viscosity)),
      gas_density_(std::move(gas_density)),
      gas_viscosity_(std::move(gas_viscosity)),
      intrinsic_permeability_models_(std::move(intrinsic_permeability_models)),
      porosity_models_(std::move(porosity_models)),
      storage_models_(std::move(storage_models)),
      capillary_pressure_models_(std::move(capillary_pressure_models)),
      relative_permeability_models_(std::move(relative_permeability_models))
{
    DBUG("Create material properties for Two-Phase flow with PP model.");
}

int TwoPhaseFlowWithPPMaterialProperties::getMaterialID(
    const std::size_t element_id) const
{
    if (material_ids_.empty())
    {
        return 0;
    }

    assert(element_id < material_ids_.size());
    return material_ids_[element_id];
}

double TwoPhaseFlowWithPPMaterialProperties::getLiquidDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return liquid_density_->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasDensity(const double p,
                                                           const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return gas_density_->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getLiquidViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return liquid_viscosity_->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return gas_viscosity_->getValue(vars);
}

Eigen::MatrixXd TwoPhaseFlowWithPPMaterialProperties::getPermeability(
    const int material_id, const double t,
    const ParameterLib::SpatialPosition& pos, const int /*dim*/) const
{
    return intrinsic_permeability_models_[material_id]->getValue(t, pos, 0, 0);
}

double TwoPhaseFlowWithPPMaterialProperties::getPorosity(
    const int material_id, const double t,
    const ParameterLib::SpatialPosition& pos, const double /*p*/,
    const double T, const double porosity_variable) const
{
    return porosity_models_[material_id]->getValue(t, pos, porosity_variable,
                                                   T);
}

double TwoPhaseFlowWithPPMaterialProperties::getSaturation(
    const int material_id, const double /*t*/,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double pc) const
{
    return capillary_pressure_models_[material_id]->getSaturation(pc);
}

double TwoPhaseFlowWithPPMaterialProperties::getCapillaryPressure(
    const int material_id, const double /*t*/,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double saturation) const
{
    return capillary_pressure_models_[material_id]->getCapillaryPressure(
        saturation);
}

double TwoPhaseFlowWithPPMaterialProperties::getSaturationDerivative(
    const int material_id, const double /*t*/,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double saturation) const
{
    const double dpcdsw =
        capillary_pressure_models_[material_id]->getdPcdS(saturation);
    return 1 / dpcdsw;
}

double TwoPhaseFlowWithPPMaterialProperties::getNonwetRelativePermeability(
    const double /*t*/, const ParameterLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    if (saturation < 0.)
    {
        return 1.0;
    }
    if (saturation > 1)
    {
        return 0.0;
    }
    return boost::math::pow<3>(1 - saturation);
}

double TwoPhaseFlowWithPPMaterialProperties::getWetRelativePermeability(
    const double /*t*/, const ParameterLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    if (saturation < 0)
    {
        return 0.0;
    }
    if (saturation > 1)
    {
        return 1.0;
    }
    return boost::math::pow<3>(saturation);
}
}  // namespace TwoPhaseFlowWithPP
}  // namespace MaterialLib
