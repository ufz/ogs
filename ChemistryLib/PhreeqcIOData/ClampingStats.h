// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <string_view>

namespace ChemistryLib
{
namespace PhreeqcIOData
{
/// Statistics for clamping negative input concentrations to zero before
/// speciation. PHREEQC rejects negative concentrations, so every negative value
/// is clamped to zero unconditionally; the warning threshold (a negative value,
/// see ChemistryLib parameter \c concentration_warning_threshold) only decides
/// whether a clamp is reported. Values in
/// \f$[\text{warning\_threshold}, 0)\f$ are treated as floating-point noise and
/// clamped silently; values \f$<\text{warning\_threshold}\f$ are clamped and
/// counted as "severe" so an aggregated warning can be emitted.
///
/// The same type is used both as a per-cell leaf (built by clamp()) and as a
/// run-level accumulator (folded by operator+=). It is an aggregate (no
/// user-declared constructor, no private data) so brace-init and `= {}` reset
/// work.
struct ClampingStats
{
    std::string_view worst_component_name;  //!< aliases Component::name
    std::size_t n_cells = 0;                //!< cells with any clamping
    std::size_t n_values = 0;               //!< total values clamped
    std::size_t n_severe_cells = 0;         //!< cells with a severe clamp
    std::size_t n_severe_values = 0;        //!< total severe values clamped
    double total_clamped_amount = 0.0;
    double worst_negative_value = 0.0;

    /// Clamps one concentration: negative values are returned as 0 and the
    /// value-level counters are updated; non-negative input is returned
    /// unchanged. A value below \c warning_threshold (which is negative) is
    /// recorded as severe and may update the worst value/name. The cell-level
    /// counters (n_cells / n_severe_cells) are set by the caller after the
    /// per-component loop.
    double clamp(double c, double warning_threshold, std::string_view name);

    /// Folds another stats object into this accumulator (identity = {}).
    ClampingStats& operator+=(ClampingStats const& other);

    /// WARN for severe clamps, DBUG for silent noise clamps. No-op when empty.
    void report(double warning_threshold) const;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
