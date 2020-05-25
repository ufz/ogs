/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include "MaterialLib/PorousMedium/Permeability/Permeability.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"

namespace ParameterLib
{
class SpatialPosition;
}

namespace ProcessLib
{
namespace RichardsComponentTransport
{

class PorousMediaProperties
{
public:
    PorousMediaProperties(
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>&&
            porosity_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>&&
            intrinsic_permeability_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>&&
            specific_storage_models,
        std::vector<std::unique_ptr<
            MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
            capillary_pressure_saturation_models,
        std::vector<
            std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
            relative_permeability_models,
        std::vector<int>&& material_ids)
        : porosity_models_(std::move(porosity_models)),
          intrinsic_permeability_models_(
              std::move(intrinsic_permeability_models)),
          specific_storage_models_(std::move(specific_storage_models)),
          capillary_pressure_saturation_models_(
              std::move(capillary_pressure_saturation_models)),
          relative_permeability_models_(
              std::move(relative_permeability_models)),
          material_ids_(std::move(material_ids))
    {
    }

    PorousMediaProperties(PorousMediaProperties&& other)
        : porosity_models_(std::move(other.porosity_models_)),
          intrinsic_permeability_models_(
              std::move(other.intrinsic_permeability_models_)),
          specific_storage_models_(std::move(other.specific_storage_models_)),
          capillary_pressure_saturation_models_(
              std::move(other.capillary_pressure_saturation_models_)),
          relative_permeability_models_(
              std::move(other.relative_permeability_models_)),
          material_ids_(other.material_ids_)
    {
    }

    MaterialLib::PorousMedium::Porosity const& getPorosity(
        double t, ParameterLib::SpatialPosition const& pos) const;

    MaterialLib::PorousMedium::Permeability const& getIntrinsicPermeability(
        double t, ParameterLib::SpatialPosition const& pos) const;

    MaterialLib::PorousMedium::Storage const& getSpecificStorage(
        double t, ParameterLib::SpatialPosition const& pos) const;

    MaterialLib::PorousMedium::CapillaryPressureSaturation const&
    getCapillaryPressureSaturationModel(
        double t, ParameterLib::SpatialPosition const& pos) const;

    MaterialLib::PorousMedium::RelativePermeability const&
    getRelativePermeability(double t,
                            ParameterLib::SpatialPosition const& pos) const;

private:
    int getMaterialID(ParameterLib::SpatialPosition const& pos) const;

private:
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        porosity_models_;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>
        intrinsic_permeability_models_;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        specific_storage_models_;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        capillary_pressure_saturation_models_;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        relative_permeability_models_;
    std::vector<int> material_ids_;
};

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
