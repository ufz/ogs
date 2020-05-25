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

namespace MaterialLib
{
namespace PorousMedium
{
int PorousMediaProperties::getMaterialID(
    ParameterLib::SpatialPosition const& pos) const
{
    return material_ids_ ? (*material_ids_)[pos.getElementID().get()] : 0;
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
}  // namespace PorousMedium
}  // namespace MaterialLib
