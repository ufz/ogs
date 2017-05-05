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
namespace RichardsFlow
{
    /** This class has a collection of material properties for Richards flow model
    *  and it makes description of the properties of unsaturated porous media
    *  i.e. the fluid density and viscosity models
    *  the relative permeability models,
    *  the capillary pressure-saturation relationships.
    *  It generally provides the computation of the PDE coefficients for Richards flow.
    */

class RichardsFlowMaterialProperties
{
public:
    using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

    RichardsFlowMaterialProperties(
        boost::optional<MeshLib::PropertyVector<int> const&> const material_ids,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            fluid_density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>&&
            fluid_viscosity,
        std::vector<Eigen::MatrixXd>&&
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

    double getRelativePermeability(const double t,
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
    double getFluidDensity(const double p, const double T) const;
    double getFluidViscosity(const double p, const double T) const;

private:
    /**
    *  Material IDs must be given as mesh element properties.
    */
    boost::optional<MeshLib::PropertyVector<int> const&> const _material_ids;

    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const _fluid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> const _fluid_viscosity;

    std::vector<Eigen::MatrixXd> const _intrinsic_permeability_models;
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
