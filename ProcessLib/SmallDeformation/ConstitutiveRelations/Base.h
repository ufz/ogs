// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/StrongType.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/KelvinVector.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::SmallDeformation
{

using namespace ProcessLib::ConstitutiveRelations;
namespace KV = MathLib::KelvinVector;

template <int DisplacementDim>
using KelvinVector = KV::KelvinVectorType<DisplacementDim>;

template <int DisplacementDim>
using KelvinMatrix = KV::KelvinMatrixType<DisplacementDim>;

template <int DisplacementDim>
using GlobalDimVector = Eigen::Vector<double, DisplacementDim>;

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

}  // namespace ProcessLib::SmallDeformation
