// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
    std::vector<std::vector<double>> strain_rate_variables_per_element;
};

using CouplingScheme = std::variant<Monolithic, Staggered>;

CouplingScheme parseCouplingScheme(BaseLib::ConfigTree const& config);
}  // namespace ProcessLib::Common::HydroMechanics
