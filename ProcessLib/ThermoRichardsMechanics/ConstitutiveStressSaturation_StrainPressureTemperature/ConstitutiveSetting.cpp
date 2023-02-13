/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstitutiveSetting.h"

#include "ProcessLib/Graph/Apply.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
void ConstitutiveSetting<DisplacementDim>::eval(
    ConstitutiveModels<DisplacementDim>& models, double const t,
    double const dt, ParameterLib::SpatialPosition const& x_position,
    MaterialPropertyLib::Medium const& medium,
    TemperatureData<DisplacementDim> const& T_data,
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    KelvinVector<DisplacementDim> const& eps_arg,
    StatefulData<DisplacementDim>& state,
    StatefulDataPrev<DisplacementDim> const& prev_state,
    MaterialStateData<DisplacementDim>& mat_state,
    ConstitutiveTempData<DisplacementDim>& tmp,
    OutputData<DisplacementDim>& out,
    ConstitutiveData<DisplacementDim>& cd) const
{
    namespace G = ProcessLib::Graph;
    constexpr auto D = DisplacementDim;

    auto const aux_data = std::tuple{SpaceTimeData{x_position, t, dt},
                                     MediaData{medium}, T_data, p_cap_data};

    auto const mat_state_tuple = std::tie(mat_state);

    // TODO will eps lag one iteration behind? (since it's not updated after
    // solving the global equation system)
    std::get<StrainData<D>>(state).eps.noalias() = eps_arg;

    G::eval(std::get<BiotModel>(models), aux_data, tmp);

    G::eval(std::get<SolidMechanicsModel<D>>(models), aux_data, tmp, state,
             prev_state, mat_state_tuple, cd);

    G::eval(
        std::get<SolidCompressibilityModel<D, SolidConstitutiveRelation<D>>>(
            models),
        aux_data, tmp, cd);

    G::eval(std::get<BishopsModel>(models), aux_data, state, tmp);
    // TODO why not ordinary state tracking?
    G::eval(std::get<BishopsPrevModel>(models), aux_data, prev_state, tmp);
    G::eval(std::get<PorosityModel<D>>(models), aux_data, tmp, state,
             prev_state);

    // TODO move to local assembler?
    {
        auto const& biot_data = std::get<BiotData>(tmp);
        auto const& poro_data = std::get<PorosityData>(state);

        if (biot_data() < poro_data.phi)
        {
            OGS_FATAL(
                "ThermoRichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                biot_data(), poro_data.phi, *x_position.getElementID(),
                *x_position.getIntegrationPoint());
        }
    }

    G::eval(std::get<LiquidDensityModel<D>>(models), aux_data, out);
    G::eval(std::get<SolidDensityModel<D>>(models), aux_data, state, out);
    G::eval(std::get<GravityModel<D>>(models), state, out, tmp, cd);
    G::eval(std::get<LiquidViscosityModel<D>>(models), aux_data, out);
    G::eval(std::get<TransportPorosityModel<D>>(models), aux_data, tmp, state,
             prev_state);
    G::eval(std::get<PermeabilityModel<D>>(models), aux_data, state, out, cd,
             tmp);
    G::eval(std::get<ThermoOsmosisModel<D>>(models), aux_data, out, cd);
    G::eval(std::get<DarcyLawModel<D>>(models), aux_data, out, tmp, cd);
    G::eval(std::get<TRMHeatStorageAndFluxModel<D>>(models), aux_data, out,
             state, tmp, cd);
    G::eval(std::get<TRMVaporDiffusionModel<D>>(models), aux_data, out, state,
             tmp, cd);

    // TODO Not needed for solid mechanics (solid thermal expansion is computed
    // by the solid material model), but for fluid expansion. This duplication
    // should be avoided in the future.
    G::eval(std::get<SolidThermalExpansionModel<D>>(models), aux_data, tmp);

    G::eval(std::get<FluidThermalExpansionModel<D>>(models), aux_data, tmp,
             state, out);
    G::eval(std::get<TRMStorageModel<D>>(models), aux_data, tmp, state, out,
             prev_state, cd);
    G::eval(std::get<EqPModel<D>>(models), aux_data, state, tmp, out, cd);
    G::eval(std::get<EqTModel<D>>(models), cd);
}

template struct ConstitutiveSetting<2>;
template struct ConstitutiveSetting<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
