// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <variant>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib::Common::HydroMechanics
{
struct Monolithic
{
};

struct Staggered
{
    /// An optional input to set an algorithmic parameter of the staggered
    /// scheme. In the HydroMechanics and the ThermoHydroMechanics processes
    /// the fixed-stress split has been implemented as staggered scheme with a
    /// stabilization parameter to be set. For more details see [user guide -
    /// conventions](https://www.opengeosys.org/docs/userguide/basics/conventions/#fixed-stress-split-for-hydro-mechanical-processes).
    double fixed_stress_stabilization_parameter;

    /// Strain rate variable of the fixed-stress split, stored per element
    /// (outer index: element ID) and per integration point (inner index).
    ///
    /// The ThermoHydroMechanics process stores the volumetric strain rate
    /// with its pressure and temperature driven parts removed,
    /// \f[ \dot{\epsilon}_v - \beta_\mathrm{FS}\,\alpha'\,c_s\,\dot p
    ///     - \beta_\mathrm{FS}\,\alpha_T\,\dot T, \f]
    /// with the stabilization parameter \f$\beta_\mathrm{FS}\f$, the
    /// (ice-modified) Biot coefficient \f$\alpha'\f$, the solid skeleton
    /// compressibility \f$c_s\f$ and the volumetric solid thermal expansion
    /// coefficient \f$\alpha_T\f$. Evaluated after the mechanical solve of
    /// every coupling iteration. The mass balance equation adds back the
    /// pressure and temperature driven parts with the current iterates.
    ///
    /// The HydroMechanics process stores the same split differently: it holds
    /// only \f$\dot p\f$ of the previous coupling iteration, and the
    /// volumetric strain rate is evaluated during assembly.
    std::vector<std::vector<double>> fixed_stress_split_variables_per_element;
};

using CouplingScheme = std::variant<Monolithic, Staggered>;

CouplingScheme parseCouplingScheme(BaseLib::ConfigTree const& config);

/// The functions below access
/// Staggered::fixed_stress_split_variables_per_element and do nothing for the
/// monolithic scheme.

/// Allocates the per-element storage. The per-integration-point storage is
/// allocated by allocateFixedStressSplitVariablesOfElement(), called from the
/// local assemblers' constructors.
void allocateFixedStressSplitVariables(CouplingScheme& coupling_scheme,
                                       std::size_t n_elements);

/// Allocates the storage of the given element, initialized with zeros.
void allocateFixedStressSplitVariablesOfElement(CouplingScheme& coupling_scheme,
                                                std::size_t element_id,
                                                unsigned n_integration_points);

/// Sets the strain rate variables of the given element from \c values, and
/// returns the number of values read.
std::size_t setFixedStressSplitVariables(CouplingScheme& coupling_scheme,
                                         std::size_t element_id,
                                         double const* values);

/// Returns the strain rate variables of the given element, or an empty vector
/// for the monolithic scheme.
std::vector<double> getFixedStressSplitVariables(
    CouplingScheme const& coupling_scheme, std::size_t element_id);
}  // namespace ProcessLib::Common::HydroMechanics
