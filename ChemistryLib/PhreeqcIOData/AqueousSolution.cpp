// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AqueousSolution.h"

#include <cmath>
#include <ostream>

#include "BaseLib/Error.h"
#include "ChemistryLib/Common/ChargeBalance.h"
#include "MeshLib/PropertyVector.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void AqueousSolution::print(std::ostream& os,
                            std::size_t const chemical_system_id) const
{
    os << "temp " << temperature << "\n";

    os << "pressure " << pressure << "\n";

    double const pH_value = -std::log10(H_plus_activity[chemical_system_id]);

    switch (charge_balance)
    {
        case ChargeBalance::pH:
            os << "pH " << pH_value << " charge" << "\n";
            os << "pe " << (*pe)[chemical_system_id] << "\n";
            break;
        case ChargeBalance::pe:
            os << "pH " << pH_value << "\n";
            os << "pe " << (*pe)[chemical_system_id] << " charge" << "\n";
            break;
        case ChargeBalance::Unspecified:
            os << "pH " << pH_value << "\n";
            os << "pe " << (*pe)[chemical_system_id] << "\n";
            break;
    }

    os << "units mol/kgw\n";

    for (auto const& component : components)
    {
        os << component.name << " " << component.amount[chemical_system_id];
        component.chemical_formula.empty()
            ? os << "\n"
            : os << " as " << component.chemical_formula << "\n";
    }

    os << "\n\n";
}

ClampingStats setAqueousSolution(std::vector<double> const& concentrations,
                                 std::size_t const chemical_system_id,
                                 AqueousSolution& aqueous_solution,
                                 double const warning_threshold)
{
    ClampingStats stats;
    auto& components = aqueous_solution.components;
    for (unsigned component_id = 0; component_id < components.size();
         ++component_id)
    {
        components[component_id].amount[chemical_system_id] =
            stats.clamp(concentrations[component_id], warning_threshold,
                        components[component_id].name);
    }
    // clamp() collects value-level counters; the per-cell counters are this
    // one cell's contribution to the run-level accumulator.
    stats.n_cells = stats.n_values > 0 ? 1 : 0;
    stats.n_severe_cells = stats.n_severe_values > 0 ? 1 : 0;

    // The transport process carries the H+ activity 10^-pH as the last entry
    // of the concentrations vector for the pH "component". Unlike component
    // concentrations it is not clamped: PHREEQC assumes 1 kg of water, and
    // water autoprotolysis always yields some H+, so a non-positive H+ activity
    // is physically impossible. Clamping it to zero would make
    // pH = -log10(activity) NaN/inf, so a non-positive value here signals a
    // transport-side bug rather than floating-point noise.
    double const h_plus_activity = concentrations.back();
    if (h_plus_activity <= 0.0)
    {
        OGS_FATAL(
            "H+ activity (10^-pH) at chemical system {:d} is {:g} <= 0, which "
            "is physically impossible (PHREEQC's 1 kg-water assumption "
            "guarantees a positive H+ activity).",
            chemical_system_id, h_plus_activity);
    }
    aqueous_solution.H_plus_activity[chemical_system_id] = h_plus_activity;
    return stats;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
