/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "TwoPhaseFlowWithPPMaterialProperties.h"
#include <logog/include/logog.hpp>
#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
TwoPhaseFlowWithPPMaterialProperties::TwoPhaseFlowWithPPMaterialProperties(
    boost::optional<MeshLib::PropertyVector<int> const&> const material_ids,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        liquid_density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        liquid_viscosity,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        gas_density,
    std::unique_ptr<MaterialLib::Fluid::FluidProperty>
        gas_viscosity,
    std::vector<Eigen::MatrixXd>
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
    : _liquid_density(std::move(liquid_density)),
      _liquid_viscosity(std::move(liquid_viscosity)),
      _gas_density(std::move(gas_density)),
      _gas_viscosity(std::move(gas_viscosity)),
      _material_ids(material_ids),
      _intrinsic_permeability_models(intrinsic_permeability_models),
      _porosity_models(std::move(porosity_models)),
      _storage_models(std::move(storage_models)),
      _capillary_pressure_models(std::move(capillary_pressure_models)),
      _relative_permeability_models(std::move(relative_permeability_models))
{
    DBUG("Create material properties for Two-Phase flow with PP model.");
}

int TwoPhaseFlowWithPPMaterialProperties::getMaterialID(
    const std::size_t element_id)
{
    if (!_material_ids)
    {
        return 0;
    }

    assert(element_id < _material_ids->size());
    return (*_material_ids)[element_id];
}

double TwoPhaseFlowWithPPMaterialProperties::getLiquidDensity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _liquid_density->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasDensity(const double p,
                                                           const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_density->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasDensityDerivative(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;

    return _gas_density->getdValue(vars,
                                   MaterialLib::Fluid::PropertyVariableType::p);
}
double TwoPhaseFlowWithPPMaterialProperties::getLiquidViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _liquid_viscosity->getValue(vars);
}

double TwoPhaseFlowWithPPMaterialProperties::getGasViscosity(
    const double p, const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _gas_viscosity->getValue(vars);
}

Eigen::MatrixXd const& TwoPhaseFlowWithPPMaterialProperties::getPermeability(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const int /*dim*/) const
{
    return _intrinsic_permeability_models[material_id];
}

double TwoPhaseFlowWithPPMaterialProperties::getPorosity(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double T, const double porosity_variable) const
{
    return _porosity_models[material_id]->getValue(porosity_variable, T);
}

double TwoPhaseFlowWithPPMaterialProperties::getNonwetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _relative_permeability_models[0]->getValue(saturation);
}

double TwoPhaseFlowWithPPMaterialProperties::getWetRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _relative_permeability_models[1]->getValue(saturation);
}

double TwoPhaseFlowWithPPMaterialProperties::getSaturation(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double pc) const
{
    return _capillary_pressure_models[material_id]->getSaturation(pc);
}
double TwoPhaseFlowWithPPMaterialProperties::getSaturationDerivative(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double saturation) const
{
    const double dpcdsw =
        _capillary_pressure_models[material_id]->getdPcdS(saturation);
    return 1 / dpcdsw;
}
}  // end of namespace
}  // end of namespace
