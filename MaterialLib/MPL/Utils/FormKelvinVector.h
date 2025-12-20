// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

extern template MathLib::KelvinVector::KelvinVectorType<2> formKelvinVector<2>(
    MaterialPropertyLib::PropertyDataType const& values);

extern template MathLib::KelvinVector::KelvinVectorType<3> formKelvinVector<3>(
    MaterialPropertyLib::PropertyDataType const& values);

}  // namespace MaterialPropertyLib
