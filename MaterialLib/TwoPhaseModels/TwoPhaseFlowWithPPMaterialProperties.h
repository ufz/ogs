/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

// TODO
// The matereial properties for two phase flow process need to be restructured
// and moved to a better place.
namespace ProcessLib
{
class SpatialPosition;
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
/** This class has a collection of material properties for two-phase flow with
* PP model
*  and it makes description of the material properties for two-phase condition,
*  i.e. the gas/liquid density and viscosity models, respectively,
*  the relative permeability models with respect to two phases,
*  the capillary pressure-saturation relationships.
*  It generally provides the computation of the PDE coefficients for two-phase
* flow.
*/
class TwoPhaseFlowWithPPMaterialProperties
{
public:
    using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

    TwoPhaseFlowWithPPMaterialProperties(
        MeshLib::PropertyVector<int> const& material_ids,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            liquid_density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            liquid_viscosity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            gas_density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            gas_viscosity,
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
            relative_permeability_models);

    int getMaterialID(const std::size_t element_id) const;

    Eigen::MatrixXd const& getPermeability(
        const int material_id,
        const double t,
        const ProcessLib::SpatialPosition& pos,
        const int dim) const;

    double getPorosity(const int material_id, const double t,
                       const ProcessLib::SpatialPosition& pos, const double p,
                       const double T, const double porosity_variable) const;
    double getSaturation(const int material_id, const double t,
                         const ProcessLib::SpatialPosition& pos, const double p,
                         const double T, const double pc) const;
    double getCapillaryPressure(const int material_id, const double t,
                                const ProcessLib::SpatialPosition& pos,
                                const double p, const double T,
                                const double saturation) const;
    double getSaturationDerivative(const int material_id, const double t,
                                   const ProcessLib::SpatialPosition& pos,
                                   const double p, const double T,
                                   const double saturation) const;
    double getLiquidDensity(const double p, const double T) const;
    double getGasDensity(const double p, const double T) const;
    double getGasViscosity(const double p, const double T) const;
    double getLiquidViscosity(const double p, const double T) const;
    double getGasDensityDerivative(double const p, double const T) const;
    double getNonwetRelativePermeability(const double t,
                                         const ProcessLib::SpatialPosition& pos,
                                         const double p, const double T,
                                         const double saturation) const;
    double getWetRelativePermeability(const double t,
                                      const ProcessLib::SpatialPosition& pos,
                                      const double p, const double T,
                                      const double saturation) const;

private:
    /** Use two phase models for different material zones.
    *  Material IDs must be given as mesh element properties.
    */
    MeshLib::PropertyVector<int> const& _material_ids;

    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const _liquid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const _liquid_viscosity;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const _gas_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const _gas_viscosity;

    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>> const
        _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>> const
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>> const
        _storage_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>> const
        _capillary_pressure_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>> const
        _relative_permeability_models;
};

}  // end of namespace
}  // end of namespace
