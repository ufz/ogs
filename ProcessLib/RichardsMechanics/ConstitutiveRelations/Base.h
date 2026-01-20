// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/StrongType.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/SpatialPosition.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::RichardsMechanics
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
        : medium{medium},
          liquid{medium.phase("AqueousLiquid")},
          solid{medium.phase("Solid")}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& liquid;
    MaterialPropertyLib::Phase const& solid;
};

using TemperatureData = BaseLib::StrongType<double, struct TemperatureDataTag>;

template <int DisplacementDim>
struct CapillaryPressureData
{
    double p_cap;
    double p_cap_prev;
    Eigen::Vector<double, DisplacementDim> grad_p_cap;
};
}  // namespace ProcessLib::RichardsMechanics
