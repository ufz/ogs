/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PorousMediaProperties.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
int PorousMediaProperties::getMaterialID(SpatialPosition const& pos) const
{
    int const element_id = pos.getElementID().get();
    return _material_ids[element_id];
}

MaterialLib::PorousMedium::Porosity const& PorousMediaProperties::getPorosity(
    double /*t*/, SpatialPosition const& pos) const
{
    return *_porosity_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Permeability const&
PorousMediaProperties::getIntrinsicPermeability(
    double /*t*/, SpatialPosition const& pos) const
{
    return *_intrinsic_permeability_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Storage const&
PorousMediaProperties::getSpecificStorage(double /*t*/,
                                          SpatialPosition const& pos) const
{
    return *_specific_storage_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::CapillaryPressureSaturation const&
PorousMediaProperties::getCapillaryPressureSaturationModel(
    double /*t*/, SpatialPosition const& pos) const
{
    return *_capillary_pressure_saturation_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::RelativePermeability const&
PorousMediaProperties::getRelativePermeability(
    double /*t*/, SpatialPosition const& pos) const
{
    return *_relative_permeability_models[getMaterialID(pos)];
}
}
}
