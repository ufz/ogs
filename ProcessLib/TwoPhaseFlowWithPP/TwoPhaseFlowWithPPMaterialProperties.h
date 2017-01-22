/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#pragma once

#include <memory>
#include <vector>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
class SpatialPosition;
namespace TwoPhaseFlowWithPP
{
    /** This class has a collection of material properties for two-phase flow with PP model
    *  and it makes description of the material properties for two-phase condition,
    *  i.e. the gas/liquid density and viscosity models, respectively,
    *  the relative permeability models with respect to two phases,
    *  the capillary pressure-saturation relationships.
    *  It generally provides the computation of the PDE coefficients for two-phase flow.
    */

class TwoPhaseFlowWithPPMaterialProperties
{
public:
    using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

    TwoPhaseFlowWithPPMaterialProperties(
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
            relative_permeability_models);

    int getMaterialID(const std::size_t element_id);

    Eigen::MatrixXd const& getPermeability(
        const int material_id,
        const double t,
        const ProcessLib::SpatialPosition& pos,
        const int dim) const;

    double getPorosity(const int material_id, const double t,
                       const ProcessLib::SpatialPosition& pos, const double p,
                       const double T, const double porosity_variable) const;

    double getNonwetRelativePermeability(const double t,
                                         const ProcessLib::SpatialPosition& pos,
                                         const double p, const double T,
                                         const double saturation) const;
    double getWetRelativePermeability(const double t,
                                      const ProcessLib::SpatialPosition& pos,
                                      const double p, const double T,
                                      const double saturation) const;
    double getSaturation(const int material_id, const double t,
                         const ProcessLib::SpatialPosition& pos, const double p,
                         const double T, const double pc) const;
    double getSaturationDerivative(const int material_id, const double t,
                                   const ProcessLib::SpatialPosition& pos,
                                   const double p, const double T,
                                   const double saturation) const;
    double getLiquidDensity(const double p, const double T) const;
    double getGasDensity(const double p, const double T) const;
    double getGasViscosity(const double p, const double T) const;
    double getLiquidViscosity(const double p, const double T) const;
    double getGasDensityDerivative(double const p, double const T) const;

protected:
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _liquid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _liquid_viscosity;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_viscosity;

    /** Use two phase models for different material zones.
    *  Material IDs must be given as mesh element properties.
    */
    boost::optional<MeshLib::PropertyVector<int> const&> const _material_ids;

    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        _capillary_pressure_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        _relative_permeability_models;
};

}  // end of namespace
}  // end of namespace
