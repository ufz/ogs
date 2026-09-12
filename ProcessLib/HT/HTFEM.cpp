// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "HTFEM.h"

#include <cmath>

#include "BaseLib/Error.h"

namespace ProcessLib
{
namespace HT
{
namespace
{
/// Checks the requirement \f$\alpha_B=1 \Rightarrow S_s=0\f$ on already
/// evaluated values and ends the run with \c OGS_FATAL if it is violated.
///
/// A specific storage that is not a number is not a violation. Where the
/// check runs before the primary variables are known, a storage property
/// depending on them evaluates to NaN, and NaN compares unequal to zero; such
/// a setting is left to the checks during assembly, where the storage has a
/// value. The same reasoning covers a Biot coefficient of NaN, which compares
/// unequal to one.
void checkBiotStorageRelationValues(ParameterLib::SpatialPosition const& pos,
                                    double const biot_coefficient,
                                    double const specific_storage)
{
    if (biot_coefficient != 1.0 || specific_storage == 0.0 ||
        std::isnan(specific_storage))
    {
        return;
    }

    OGS_FATAL(
        "At {} the Biot coefficient evaluates to 1.0, which requires the "
        "specific storage of the solid phase to be 0.0, but it evaluates to "
        "{:g}.",
        pos, specific_storage);
}
}  // namespace

double evalEffectiveThermalExpansivity(
    double const t, double const dt, ParameterLib::SpatialPosition const& pos,
    MaterialPropertyLib::VariableArray const& vars,
    MaterialPropertyLib::Medium const& medium,
    MaterialPropertyLib::Phase const& liquid_phase,
    MaterialPropertyLib::Phase const& solid_phase,
    bool const has_solid_thermal_expansivity, double const specific_storage)
{
    double const dfluid_density_dT =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                vars, MaterialPropertyLib::Variable::temperature, pos, t, dt);

    double const fluid_density = vars.density;
    double const porosity = vars.porosity;
    double const fluid_thermal_expansivity =
        -porosity * dfluid_density_dT / fluid_density;

    if (!has_solid_thermal_expansivity)
    {
        return fluid_thermal_expansivity;
    }

    double const linear_solid_thermal_expansivity =
        solid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_expansivity)
            .template value<double>(vars, pos, t, dt);
    double const biot_coefficient =
        medium.property(MaterialPropertyLib::PropertyType::biot_coefficient)
            .template value<double>(vars, pos, t, dt);

    // Both properties may vary in space and time, so the evaluated values are
    // the only ones that can be compared. The state independent cases are
    // caught earlier by checkBiotStorageRelation().
    checkBiotStorageRelationValues(pos, biot_coefficient, specific_storage);

    return fluid_thermal_expansivity + 3.0 * (biot_coefficient - porosity) *
                                           linear_solid_thermal_expansivity;
}

void checkBiotStorageRelation(double const t, double const dt,
                              ParameterLib::SpatialPosition const& pos,
                              MaterialPropertyLib::VariableArray const& vars,
                              MaterialPropertyLib::Medium const& medium,
                              MaterialPropertyLib::Phase const& solid_phase)
{
    double const biot_coefficient =
        medium.property(MaterialPropertyLib::PropertyType::biot_coefficient)
            .template value<double>(vars, pos, t, dt);
    double const specific_storage =
        solid_phase.property(MaterialPropertyLib::PropertyType::storage)
            .template value<double>(vars, pos, t, dt);

    checkBiotStorageRelationValues(pos, biot_coefficient, specific_storage);
}
}  // namespace HT
}  // namespace ProcessLib
