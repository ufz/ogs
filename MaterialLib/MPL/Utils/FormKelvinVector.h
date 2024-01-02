/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 12, 2021, 4:34 PM
 */

#pragma once

#include <Eigen/Core>
#include <variant>

#include "MaterialLib/MPL/Property.h"
#include "MathLib/KelvinVector.h"

namespace MaterialPropertyLib
{
/**
 * \brief A function to form a Kelvin vector from strain or stress alike
 * property like thermal expansivity for thermal strain.
 *
 * It takes either a scalar number for isotropic thermal expansion or a
 * three element vector or a 3 x 3 matrix for anisotropic properties, to
 * get a Kelvin vector for strain or stress.
 *
 * \param values: e.g., Thermal expansivity, which can be scalar number, a three
 *                element vector or a 3 x 3 matrix.
 * \return        Kelvin vector type property.
 */
template <int GlobalDim>
MathLib::KelvinVector::KelvinVectorType<GlobalDim> formKelvinVector(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
