// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ClampingStats.h"

#include "BaseLib/Logging.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
double ClampingStats::clamp(double const c, double const warning_threshold,
                            std::string_view const name)
{
    if (c >= 0.0)
    {
        return c;
    }
    if (c < warning_threshold)
    {
        ++n_severe_values;
        if (c < worst_negative_value)
        {
            worst_negative_value = c;
            worst_component_name = name;
        }
    }
    ++n_values;
    total_clamped_amount += -c;
    return 0.0;
}

ClampingStats& ClampingStats::operator+=(ClampingStats const& other)
{
    n_cells += other.n_cells;
    n_values += other.n_values;
    n_severe_cells += other.n_severe_cells;
    n_severe_values += other.n_severe_values;
    total_clamped_amount += other.total_clamped_amount;
    if (other.worst_negative_value < worst_negative_value)
    {
        worst_negative_value = other.worst_negative_value;
        worst_component_name = other.worst_component_name;
    }
    return *this;
}

void ClampingStats::report(double const warning_threshold) const
{
    if (n_severe_values > 0)
    {
        WARN(
            "Clamped {:d} concentration value(s) more negative than the "
            "warning threshold ({:g} mol/kgw) to zero across {:d} chemical "
            "system(s) before speciation; worst: '{:s}' = {:g} mol/kgw.",
            n_severe_values, warning_threshold, n_severe_cells,
            worst_component_name, worst_negative_value);
    }
    if (n_values > 0)
    {
        DBUG(
            "PhreeqcIO: clamped {:d} negative concentration value(s) to zero "
            "across {:d} chemical system(s) before speciation (total clamped: "
            "{:g} mol/kgw).",
            n_values, n_cells, total_clamped_amount);
    }
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
