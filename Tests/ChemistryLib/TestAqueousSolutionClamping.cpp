// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include "ChemistryLib/Common/ChargeBalance.h"
#include "ChemistryLib/PhreeqcIOData/AqueousSolution.h"
#include "ChemistryLib/PhreeqcIOData/ClampingStats.h"

using namespace ChemistryLib;
using namespace ChemistryLib::PhreeqcIOData;

namespace
{
// Builds an AqueousSolution with the given component names. Each component's
// amount vector and the H+ activity vector are sized for num_systems chemical
// systems so that setAqueousSolution() can index them by chemical_system_id.
AqueousSolution makeAqueousSolution(std::vector<std::string> const& names,
                                    std::size_t const num_systems)
{
    std::vector<Component> components;
    for (auto const& name : names)
    {
        Component component(name, /*chemical_formula=*/"");
        component.amount.assign(num_systems, -999.0);  // sentinel
        components.push_back(std::move(component));
    }

    AqueousSolution aqueous_solution(
        /*fixing_pe=*/false, /*temperature=*/298.15, /*pressure=*/1.0e5,
        /*pe=*/nullptr, /*pe0=*/4.0, std::move(components),
        ChargeBalance::Unspecified);
    aqueous_solution.H_plus_activity.assign(num_systems, -999.0);
    return aqueous_solution;
}

// Signed warning threshold: clamps below this value are reported as severe,
// values in [warning_threshold, 0) are clamped silently.
constexpr double warning_threshold = -1e-12;
}  // namespace

// All concentrations non-negative: nothing is clamped, amounts are written
// through unchanged and the H+ activity is taken from the last entry.
TEST(ChemistryLibAqueousSolutionClamping, NoNegativesLeavesValuesUnchanged)
{
    auto aqueous_solution = makeAqueousSolution({"Ca", "Cl", "H"}, 1);
    std::vector<double> const concentrations = {0.1, 0.2, 1e-7};

    auto const stats = setAqueousSolution(concentrations, 0, aqueous_solution,
                                          warning_threshold);

    EXPECT_EQ(0u, stats.n_values);
    EXPECT_EQ(0u, stats.n_severe_values);
    EXPECT_EQ(0u, stats.n_cells);
    EXPECT_EQ(0u, stats.n_severe_cells);
    EXPECT_EQ(0.0, stats.total_clamped_amount);
    EXPECT_TRUE(stats.worst_component_name.empty());

    EXPECT_EQ(0.1, aqueous_solution.components[0].amount[0]);
    EXPECT_EQ(0.2, aqueous_solution.components[1].amount[0]);
    EXPECT_EQ(1e-7, aqueous_solution.components[2].amount[0]);
    EXPECT_EQ(1e-7, aqueous_solution.H_plus_activity[0]);
}

// A negative within [warning_threshold, 0) is treated as floating-point noise:
// it is clamped to zero and counted, but not flagged as a severe clamp.
TEST(ChemistryLibAqueousSolutionClamping, TinyNegativeIsClampedSilently)
{
    auto aqueous_solution = makeAqueousSolution({"Ca", "H"}, 1);
    std::vector<double> const concentrations = {-1e-15, 1e-7};

    auto const stats = setAqueousSolution(concentrations, 0, aqueous_solution,
                                          warning_threshold);

    EXPECT_EQ(1u, stats.n_values);
    EXPECT_EQ(0u, stats.n_severe_values);
    EXPECT_EQ(1u, stats.n_cells);
    EXPECT_EQ(0u, stats.n_severe_cells);
    EXPECT_DOUBLE_EQ(1e-15, stats.total_clamped_amount);
    EXPECT_TRUE(stats.worst_component_name.empty());

    EXPECT_EQ(0.0, aqueous_solution.components[0].amount[0]);
}

// A negative below warning_threshold is clamped to zero and reported as a
// severe clamp, recording the offending component and its value.
TEST(ChemistryLibAqueousSolutionClamping, SevereNegativeIsReported)
{
    auto aqueous_solution = makeAqueousSolution({"Ca", "H"}, 1);
    std::vector<double> const concentrations = {-0.05, 1e-7};

    auto const stats = setAqueousSolution(concentrations, 0, aqueous_solution,
                                          warning_threshold);

    EXPECT_EQ(1u, stats.n_values);
    EXPECT_EQ(1u, stats.n_severe_values);
    EXPECT_EQ(1u, stats.n_cells);
    EXPECT_EQ(1u, stats.n_severe_cells);
    EXPECT_DOUBLE_EQ(0.05, stats.total_clamped_amount);
    EXPECT_DOUBLE_EQ(-0.05, stats.worst_negative_value);
    EXPECT_EQ("Ca", stats.worst_component_name);

    EXPECT_EQ(0.0, aqueous_solution.components[0].amount[0]);
}

// clamp() returns non-negative input unchanged and records no clamping.
TEST(ChemistryLibClampingStats, ClampPassesNonNegativeThrough)
{
    ClampingStats stats;
    EXPECT_EQ(0.5, stats.clamp(0.5, warning_threshold, "Ca"));
    EXPECT_EQ(0.0, stats.clamp(0.0, warning_threshold, "Ca"));

    EXPECT_EQ(0u, stats.n_values);
    EXPECT_EQ(0u, stats.n_severe_values);
    EXPECT_EQ(0.0, stats.total_clamped_amount);
    EXPECT_TRUE(stats.worst_component_name.empty());
}

// clamp() of a tiny negative returns zero, counts it, but leaves it non-severe
// (no worst value/name recorded).
TEST(ChemistryLibClampingStats, ClampTinyNegativeIsNotSevere)
{
    ClampingStats stats;
    EXPECT_EQ(0.0, stats.clamp(-1e-15, warning_threshold, "Ca"));

    EXPECT_EQ(1u, stats.n_values);
    EXPECT_EQ(0u, stats.n_severe_values);
    EXPECT_DOUBLE_EQ(1e-15, stats.total_clamped_amount);
    EXPECT_EQ(0.0, stats.worst_negative_value);
    EXPECT_TRUE(stats.worst_component_name.empty());
}

// clamp() of a value below the threshold returns zero and records it as severe
// together with the worst value and component name.
TEST(ChemistryLibClampingStats, ClampSevereNegativeRecordsWorst)
{
    ClampingStats stats;
    EXPECT_EQ(0.0, stats.clamp(-0.05, warning_threshold, "Ca"));

    EXPECT_EQ(1u, stats.n_values);
    EXPECT_EQ(1u, stats.n_severe_values);
    EXPECT_DOUBLE_EQ(0.05, stats.total_clamped_amount);
    EXPECT_DOUBLE_EQ(-0.05, stats.worst_negative_value);
    EXPECT_EQ("Ca", stats.worst_component_name);
}

// operator+= sums all counters and keeps the most negative worst value along
// with the component name that produced it.
TEST(ChemistryLibClampingStats, AccumulateKeepsMostNegativeWorst)
{
    ClampingStats a;
    a.clamp(-0.05, warning_threshold, "Ca");
    a.n_cells = 1;
    a.n_severe_cells = 1;

    ClampingStats b;
    b.clamp(-0.10, warning_threshold, "Cl");
    b.n_cells = 1;
    b.n_severe_cells = 1;

    a += b;

    EXPECT_EQ(2u, a.n_values);
    EXPECT_EQ(2u, a.n_severe_values);
    EXPECT_EQ(2u, a.n_cells);
    EXPECT_EQ(2u, a.n_severe_cells);
    EXPECT_DOUBLE_EQ(0.15, a.total_clamped_amount);
    EXPECT_DOUBLE_EQ(-0.10, a.worst_negative_value);
    EXPECT_EQ("Cl", a.worst_component_name);
}

// A non-severe leaf (worst value 0, empty name) must not overwrite an existing
// worst value/name when folded in.
TEST(ChemistryLibClampingStats, AccumulateNonSevereDoesNotOverwriteWorst)
{
    ClampingStats a;
    a.clamp(-0.05, warning_threshold, "Ca");
    a.n_cells = 1;
    a.n_severe_cells = 1;

    ClampingStats b;
    b.clamp(-1e-15, warning_threshold, "Cl");  // not severe
    b.n_cells = 1;

    a += b;

    EXPECT_EQ(2u, a.n_values);
    EXPECT_EQ(1u, a.n_severe_values);
    EXPECT_DOUBLE_EQ(-0.05, a.worst_negative_value);
    EXPECT_EQ("Ca", a.worst_component_name);
}

// report() on an empty accumulator is a no-op and must not crash.
TEST(ChemistryLibClampingStats, ReportOnEmptyDoesNotCrash)
{
    ClampingStats{}.report(warning_threshold);
    SUCCEED();
}
