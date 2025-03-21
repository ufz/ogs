/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Porosity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void PorosityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    SolidCompressibilityData const& solid_compressibility_data,
    SaturationData const& S_L_data,
    PrevState<SaturationData> const& S_L_prev_data,
    BishopsData const& bishops_data,
    PrevState<BishopsData> const& bishops_data_prev,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    StrainData<DisplacementDim> const& eps_data,
    PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
    PrevState<PorosityData> const& poro_prev_data, PorosityData& out) const
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;
    // TODO Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables.grain_compressibility = solid_compressibility_data.beta_SR;

    variables.liquid_saturation = S_L_data.S_L;
    variables_prev.liquid_saturation = S_L_prev_data->S_L;

    variables.effective_pore_pressure =
        -bishops_data.chi_S_L * p_cap_data.p_cap;

    // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    variables_prev.effective_pore_pressure =
        -bishops_data_prev->chi_S_L * p_cap_data.p_cap_prev;

    // Used in MaterialLib/MPL/Properties/PorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/TransportPorosityFromMassBalance.cpp
    // and MaterialLib/MPL/Properties/StrainDependentPermeability.cpp
    // Set volumetric strain rate for the general case without swelling.
    variables.volumetric_strain = Invariants::trace(eps_data.eps);
    variables_prev.volumetric_strain = Invariants::trace(eps_prev_data->eps);

    // Porosity update
    variables_prev.porosity = poro_prev_data->phi;
    out.phi = media_data.medium.property(MPL::PropertyType::porosity)
                  .template value<double>(variables, variables_prev, x_t.x,
                                          x_t.t, x_t.dt);
}

template struct PorosityModel<2>;
template struct PorosityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
