/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"

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
        std::vector<std::unique_ptr<
            MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
            capillary_pressure_saturation_models,
        std::vector<int>&& material_ids)
        : _capillary_pressure_saturation_models(
              std::move(capillary_pressure_saturation_models)),
          _material_ids(std::move(material_ids))
    {
    }

    PorousMediaProperties(PorousMediaProperties&& other)
        : _capillary_pressure_saturation_models(
              std::move(other._capillary_pressure_saturation_models)),
          _material_ids(other._material_ids)
    {
    }

    MaterialLib::PorousMedium::CapillaryPressureSaturation const&
    getCapillaryPressureSaturationModel(
        double t, ParameterLib::SpatialPosition const& pos) const;

private:
    int getMaterialID(ParameterLib::SpatialPosition const& pos) const;

private:
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        _capillary_pressure_saturation_models;
    std::vector<int> _material_ids;
};

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
