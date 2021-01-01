/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    return _material_ids ? (*_material_ids)[pos.getElementID().get()] : 0;
}

MaterialLib::PorousMedium::Porosity const& PorousMediaProperties::getPorosity(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *_porosity_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Permeability const&
PorousMediaProperties::getIntrinsicPermeability(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *_intrinsic_permeability_models[getMaterialID(pos)];
}

MaterialLib::PorousMedium::Storage const&
PorousMediaProperties::getSpecificStorage(
    double /*t*/, ParameterLib::SpatialPosition const& pos) const
{
    return *_specific_storage_models[getMaterialID(pos)];
}
}  // namespace PorousMedium
}  // namespace MaterialLib
