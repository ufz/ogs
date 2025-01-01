/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
namespace ConstitutiveRelations
{
using namespace ProcessLib::ConstitutiveRelations;
namespace KV = MathLib::KelvinVector;

template <int DisplacementDim>
using KelvinVector = KV::KelvinVectorType<DisplacementDim>;

template <int DisplacementDim>
using KelvinMatrix = KV::KelvinMatrixType<DisplacementDim>;

template <int DisplacementDim>
using GlobalDimMatrix =
    Eigen::Matrix<double, DisplacementDim, DisplacementDim, Eigen::RowMajor>;

template <int DisplacementDim>
using GlobalDimVector = Eigen::Vector<double, DisplacementDim>;

struct MediaData
{
    explicit MediaData(MaterialPropertyLib::Medium const& medium)
        : medium{medium}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& solid = medium.phase("Solid");
    MaterialPropertyLib::Phase const& liquid = medium.phase("AqueousLiquid");
    MaterialPropertyLib::Phase const& gas = medium.phase("Gas");
};

struct TemperatureData
{
    double T = nan;
    double T_prev = nan;
};

using ReferenceTemperatureData =
    BaseLib::StrongType<double, struct ReferenceTemperatureTag>;
using GasPressureData = BaseLib::StrongType<double, struct GasPressureTag>;
using CapillaryPressureData =
    BaseLib::StrongType<double, struct CapillaryPressureTag>;
template <int DisplacementDim>
using GasPressureGradientData =
    BaseLib::StrongType<GlobalDimVector<DisplacementDim>,
                        struct GasPressureGradientTag>;
template <int DisplacementDim>
using CapillaryPressureGradientData =
    BaseLib::StrongType<GlobalDimVector<DisplacementDim>,
                        struct CapillaryPressureGradientTag>;
template <int DisplacementDim>
using TemperatureGradientData =
    BaseLib::StrongType<GlobalDimVector<DisplacementDim>,
                        struct TemperatureGradientTag>;

template <int DisplacementDim>
using SpecificBodyForceData =
    BaseLib::StrongType<GlobalDimVector<DisplacementDim>,
                        struct SpecificBodyForceTag>;

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
