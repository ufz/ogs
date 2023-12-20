/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BaseLib/StrongType.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"

namespace ProcessLib::TH2M
{
using namespace ProcessLib::ConstitutiveRelations;

template <int DisplacementDim>
using KelvinVector = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

template <int DisplacementDim>
using KelvinMatrix = MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

template <int DisplacementDim>
using GlobalDimMatrix =
    Eigen::Matrix<double, DisplacementDim, DisplacementDim, Eigen::RowMajor>;

/// Convenience alias for not a number.
static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

/// Used to set a Kelvin matrix to all not-a-number.
template <int DisplacementDim>
constexpr KelvinMatrix<DisplacementDim> KMnan()
{
    return KelvinMatrix<DisplacementDim>::Constant(nan);
}

/// Used to set a Kelvin vector to all zero.
template <int DisplacementDim>
constexpr KelvinVector<DisplacementDim> KVzero()
{
    return KelvinVector<DisplacementDim>::Zero();
}

struct MediaData
{
    explicit MediaData(MaterialPropertyLib::Medium const& medium)
        : medium{medium}, solid{medium.phase("Solid")}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& solid;
};

struct TemperatureData
{
    double T;
    double T_prev;
};

using CapillaryPressureData =
    BaseLib::StrongType<double, struct CapillaryPressureTag>;

}  // namespace ProcessLib::TH2M
