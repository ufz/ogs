/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BaseLib/StrongType.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/Tensor.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::LargeDeformation
{
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

/// Represents a previous state of type T.
template <typename T>
struct PrevState
{
    PrevState() = default;
    explicit PrevState(T const& t) : t{t} {}
    explicit PrevState(T&& t) : t{std::move(t)} {}

    PrevState<T>& operator=(T const& u)
    {
        t = u;
        return *this;
    }

    PrevState<T>& operator=(T&& u)
    {
        t = std::move(u);
        return *this;
    }

    T& operator*() { return t; }
    T const& operator*() const { return t; }

    T* operator->() { return &t; }
    T const* operator->() const { return &t; }

private:
    T t;
};

struct SpaceTimeData
{
    ParameterLib::SpatialPosition x;
    double t;
    double dt;
};

struct MediaData
{
    explicit MediaData(MaterialPropertyLib::Medium const& medium)
        : medium{medium}, solid{medium.phase("Solid")}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& solid;
};

using Temperature = BaseLib::StrongType<double, struct TemperatureTag>;

template <int DisplacementDim>
struct DeformationGradientData
{
    // TODO Move initialization to the local assembler.
    MaterialPropertyLib::Tensor<DisplacementDim> deformation_gradient =
        MaterialPropertyLib::Tensor<DisplacementDim>::Zero();

    static auto reflect()
    {
        using Self = DeformationGradientData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName(
            "deformation_gradient", &Self::deformation_gradient);
    }
};

template <int DisplacementDim>
struct StrainData
{
    // TODO Move initialization to the local assembler.
    KelvinVector<DisplacementDim> eps = KVnan<DisplacementDim>();

    static auto reflect()
    {
        using Self = StrainData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("epsilon", &Self::eps);
    }
};

}  // namespace ProcessLib::LargeDeformation
