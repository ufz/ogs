/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{

using namespace ProcessLib::ConstitutiveRelations;

template <int DisplacementDim>
using KelvinVector = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

template <int DisplacementDim>
using KelvinMatrix = MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

template <int DisplacementDim>
using GlobalDimVector = Eigen::Vector<double, DisplacementDim>;

template <int DisplacementDim>
using GlobalDimMatrix =
    Eigen::Matrix<double, DisplacementDim, DisplacementDim, Eigen::RowMajor>;

/// Convenience alias for not a number.
static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

/// Used to set a Kelvin vector to all not-a-number.
template <int DisplacementDim>
constexpr KelvinVector<DisplacementDim> KVnan()
{
    return KelvinVector<DisplacementDim>::Constant(nan);
}

/// Used to set a Kelvin matrix to all not-a-number.
template <int DisplacementDim>
constexpr KelvinMatrix<DisplacementDim> KMnan()
{
    return KelvinMatrix<DisplacementDim>::Constant(nan);
}

/// Used to set a D dimensional vector to all not-a-number.
template <int D>
constexpr GlobalDimVector<D> DVnan()
{
    return GlobalDimVector<D>::Constant(nan);
}

/// Used to set a D x D matrix to all not-a-number.
template <int D>
constexpr GlobalDimMatrix<D> DMnan()
{
    return GlobalDimMatrix<D>::Constant(nan);
}

/// Used to set a Kelvin vector to all zero.
template <int DisplacementDim>
constexpr KelvinVector<DisplacementDim> KVzero()
{
    return KelvinVector<DisplacementDim>::Zero();
}

/// Used to set a Kelvin matrix to all zero.
template <int DisplacementDim>
constexpr KelvinMatrix<DisplacementDim> KMzero()
{
    return KelvinMatrix<DisplacementDim>::Zero();
}

struct MediaData
{
    explicit MediaData(MaterialPropertyLib::Medium const& medium)
        : medium{medium},
          liquid{medium.phase("AqueousLiquid")},
          solid{medium.phase("Solid")}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& liquid;
    MaterialPropertyLib::Phase const& solid;
};

template <int DisplacementDim>
struct TemperatureData
{
    double T;
    double T_prev;
    Eigen::Vector<double, DisplacementDim> grad_T;
};

template <int DisplacementDim>
struct CapillaryPressureData
{
    double p_cap;
    double p_cap_prev;
    Eigen::Vector<double, DisplacementDim> grad_p_cap;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
