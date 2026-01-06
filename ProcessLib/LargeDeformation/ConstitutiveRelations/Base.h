// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/Tensor.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::LargeDeformation
{

using namespace ProcessLib::ConstitutiveRelations;

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

struct MediaData
{
    explicit MediaData(MaterialPropertyLib::Medium const& medium)
        : medium{medium}, solid{medium.phase("Solid")}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& solid;
};

template <int DisplacementDim>
struct DeformationGradientData
{
    // TODO Move initialization to the local assembler.
    MaterialPropertyLib::Tensor<DisplacementDim> deformation_gradient =
        MaterialPropertyLib::Tensor<DisplacementDim>::Zero();
    double volume_ratio = 0;

    static auto reflect()
    {
        using Self = DeformationGradientData<DisplacementDim>;

        return std::tuple{
            ProcessLib::Reflection::makeReflectionData(
                "deformation_gradient", &Self::deformation_gradient),
            ProcessLib::Reflection::makeReflectionData("volume_ratio",
                                                       &Self::volume_ratio)};
    }
};
}  // namespace ProcessLib::LargeDeformation
