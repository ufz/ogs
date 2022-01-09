/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 12, 2021, 4:34 PM
 */

#pragma once

#include <Eigen/Dense>
#include <variant>

#include "MaterialLib/MPL/Property.h"
#include "MathLib/KelvinVector.h"

namespace MaterialPropertyLib
{
/**
 * \brief A function to form a Kelvin vector from thermal expansivity for
 * thermal strain.
 *
 * It takes the thermal expansivity, either a scalar number for isotropic
 * thermal expansion or a three element vector for anisotropic thermal
 * expansion, to get a Kelvin vector for thermal strain.
 *
 * @param values: Thermal expansivity, which can be scalar number or a three
 *                element vector.
 * @return        Thermal expansivity in Kelvin vector type.
 */
template <int GlobalDim>
MathLib::KelvinVector::KelvinVectorType<GlobalDim>
formKelvinVectorFromThermalExpansivity(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
