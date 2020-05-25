/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PorousMediaProperties.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
int PorousMediaProperties::getMaterialID(
    ParameterLib::SpatialPosition const& pos) const
{
    int const element_id = pos.getElementID().get();
    return material_ids_[element_id];
}

MaterialLib::PorousMedium::Porosity const& PorousMediaProperties::getPorosity(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *porosity_models_[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Permeability const&
PorousMediaProperties::getIntrinsicPermeability(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *intrinsic_permeability_models_[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Storage const&
PorousMediaProperties::getSpecificStorage(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *specific_storage_models_[getMaterialID(pos)];
}

MaterialLib::PorousMedium::CapillaryPressureSaturation const&
PorousMediaProperties::getCapillaryPressureSaturationModel(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *capillary_pressure_saturation_models_[getMaterialID(pos)];
}

MaterialLib::PorousMedium::RelativePermeability const&
PorousMediaProperties::getRelativePermeability(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *relative_permeability_models_[getMaterialID(pos)];
}
}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
