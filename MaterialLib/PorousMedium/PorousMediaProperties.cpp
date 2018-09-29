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

namespace MaterialLib
{
namespace PorousMedium
{
int PorousMediaProperties::getMaterialID(
    ProcessLib::SpatialPosition const& pos) const
{
    int const element_id = pos.getElementID().get();
    return _material_ids[element_id];
}

MaterialLib::PorousMedium::Porosity const& PorousMediaProperties::getPorosity(
    double /*t*/, ProcessLib::SpatialPosition const& pos) const
{
    return *_porosity_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Permeability const&
PorousMediaProperties::getIntrinsicPermeability(
    double /*t*/, ProcessLib::SpatialPosition const& pos) const
{
    return *_intrinsic_permeability_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Storage const&
PorousMediaProperties::getSpecificStorage(
    double /*t*/, ProcessLib::SpatialPosition const& pos) const
{
    return *_specific_storage_models[getMaterialID(pos)];
}
}
}
