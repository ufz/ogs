/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "RichardsFlowMaterialProperties.h"

#include <logog/include/logog.hpp>

#include "BaseLib/reorderVector.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

namespace ProcessLib
{
namespace RichardsFlow
{
RichardsFlowMaterialProperties::RichardsFlowMaterialProperties(
    boost::optional<MeshLib::PropertyVector<int> const&> const material_ids,
    std::unique_ptr<MaterialLib::Fluid::FluidProperties>&& fluid_properties,
    std::vector<Eigen::MatrixXd>&& intrinsic_permeability_models,
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
    : _material_ids(material_ids),
      _fluid_properties(std::move(fluid_properties)),
      _intrinsic_permeability_models(intrinsic_permeability_models),
      _porosity_models(std::move(porosity_models)),
      _storage_models(std::move(storage_models)),
      _capillary_pressure_models(std::move(capillary_pressure_models)),
      _relative_permeability_models(std::move(relative_permeability_models))
{
    DBUG("Create material properties for Richards flow.");
}

int RichardsFlowMaterialProperties::getMaterialID(const std::size_t element_id)
{
    if (!_material_ids)
    {
        return 0;
    }

    assert(element_id < _material_ids->size());
    return (*_material_ids)[element_id];
}

double RichardsFlowMaterialProperties::getFluidDensity(const double p,
                                                       const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _fluid_properties->getValue(
        MaterialLib::Fluid::FluidPropertyType::Density, vars);
}

double RichardsFlowMaterialProperties::getFluidViscosity(const double p,
                                                         const double T) const
{
    ArrayType vars;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = T;
    vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = p;
    return _fluid_properties->getValue(
        MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
}

Eigen::MatrixXd const& RichardsFlowMaterialProperties::getPermeability(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const int /*dim*/) const
{
    return _intrinsic_permeability_models[material_id];
}

double RichardsFlowMaterialProperties::getPorosity(
    const int material_id, const double t,
    const ProcessLib::SpatialPosition& pos, const double /*p*/,
    const double T, const double porosity_variable) const
{
    return _porosity_models[material_id]->getValue(t, pos, porosity_variable,
                                                   T);
}

double RichardsFlowMaterialProperties::getStorage(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double storage_variable) const
{
    // \todo getValue() can be extended for non
    // constant storage model
    return _storage_models[material_id]->getValue(storage_variable);
}

double RichardsFlowMaterialProperties::getRelativePermeability(
    const double /*t*/, const ProcessLib::SpatialPosition& /*pos*/,
    const double /*p*/, const double /*T*/, const double saturation) const
{
    return _relative_permeability_models[0]->getValue(saturation);
}

double RichardsFlowMaterialProperties::getSaturation(
    const int material_id, const double /*t*/,
    const ProcessLib::SpatialPosition& /*pos*/, const double /*p*/,
    const double /*T*/, const double pc) const
{
    return _capillary_pressure_models[material_id]->getSaturation(pc);
}
double RichardsFlowMaterialProperties::getSaturationDerivative(
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
