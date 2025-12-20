// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/StrongType.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/PropertyType.h"
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
    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& solid = medium.phase("Solid");
    MaterialPropertyLib::Phase const& liquid = medium.phase("AqueousLiquid");
    MaterialPropertyLib::Phase const& gas = medium.phase("Gas");

    MaterialPropertyLib::Property const& saturation_prop =
        medium[MaterialPropertyLib::PropertyType::saturation];
    MaterialPropertyLib::Property const& permeability_prop =
        medium[MaterialPropertyLib::PropertyType::permeability];
    MaterialPropertyLib::Property const& relative_permeability_prop =
        medium[MaterialPropertyLib::PropertyType::relative_permeability];
    MaterialPropertyLib::Property const&
        relative_permeability_nonwetting_phase_prop =
            medium[MaterialPropertyLib::PropertyType::
                       relative_permeability_nonwetting_phase];
    MaterialPropertyLib::Property const& thermal_conductivity_prop =
        medium[MaterialPropertyLib::PropertyType::thermal_conductivity];
    MaterialPropertyLib::Property const& specific_heat_capacity_prop =
        solid[MaterialPropertyLib::PropertyType::specific_heat_capacity];
    MaterialPropertyLib::Property const& specific_heat_capacity_solid =
        solid[MaterialPropertyLib::PropertyType::specific_heat_capacity];
    MaterialPropertyLib::Property const* transport_porosity_prop =
        medium.hasProperty(
            MaterialPropertyLib::PropertyType::transport_porosity)
            ? &medium[MaterialPropertyLib::PropertyType::transport_porosity]
            : nullptr;
    MaterialPropertyLib::Property const& viscosity_liquid =
        liquid[MaterialPropertyLib::PropertyType::viscosity];
    MaterialPropertyLib::Property const& viscosity_gas =
        gas[MaterialPropertyLib::PropertyType::viscosity];
    MaterialPropertyLib::Property const& density_solid =
        solid[MaterialPropertyLib::PropertyType::density];
    MaterialPropertyLib::Property const& density_liquid =
        liquid[MaterialPropertyLib::PropertyType::density];
    MaterialPropertyLib::Property const& density_gas =
        gas[MaterialPropertyLib::PropertyType::density];
    MaterialPropertyLib::Property const& bishops_effective_stress_prop =
        medium[MaterialPropertyLib::PropertyType::bishops_effective_stress];
    MaterialPropertyLib::Property const& porosity_prop =
        medium[MaterialPropertyLib::PropertyType::porosity];
    MaterialPropertyLib::Property const& thermal_expansivity_prop =
        solid[MaterialPropertyLib::PropertyType::thermal_expansivity];
    MaterialPropertyLib::Property const& thermal_expansivity_solid =
        solid[MaterialPropertyLib::PropertyType::thermal_expansivity];
    MaterialPropertyLib::Property const& biot_coefficient_prop =
        medium[MaterialPropertyLib::PropertyType::biot_coefficient];
    MaterialPropertyLib::Property const* tortuosity_prop =
        medium.hasProperty(MaterialPropertyLib::PropertyType::tortuosity)
            ? &medium[MaterialPropertyLib::PropertyType::tortuosity]
            : nullptr;
    MaterialPropertyLib::Property const* swelling_stress_rate_solid =
        solid.hasProperty(
            MaterialPropertyLib::PropertyType::swelling_stress_rate)
            ? &solid[MaterialPropertyLib::PropertyType::swelling_stress_rate]
            : nullptr;
};

struct TemperatureData
{
    double T = nan;
    double T_prev = nan;
};

struct GasPressureData
{
    double pG = nan;
    double pG_prev = nan;
};

struct CapillaryPressureData
{
    double pCap = nan;
    double pCap_prev = nan;
};

using ReferenceTemperatureData =
    BaseLib::StrongType<double, struct ReferenceTemperatureTag>;

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
