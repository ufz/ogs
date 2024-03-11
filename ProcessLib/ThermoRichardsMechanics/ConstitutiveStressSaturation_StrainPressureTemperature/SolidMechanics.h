/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/EquivalentPlasticStrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/MaterialState.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Saturation.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidMechanicsDataStateless.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TotalStressData.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
struct SolidMechanicsModel
{
    explicit SolidMechanicsModel(
        SolidConstitutiveRelation<DisplacementDim> const& solid_material)
        : solid_material_(solid_material),
          tangent_operator_blocks_view_{
              solid_material.template createTangentOperatorBlocksView()}
    {
    }

    void eval(const SpaceTimeData& x_t,
              TemperatureData<DisplacementDim> const& T_data,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              StrainData<DisplacementDim> const& eps_data,
              PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
              MaterialStateData<DisplacementDim>& mat_state,
              PrevState<TotalStressData<DisplacementDim>> const&
                  total_stress_data_prev,
              TotalStressData<DisplacementDim>& total_stress_data,
              EquivalentPlasticStrainData& equiv_plast_strain_data,
              SolidMechanicsDataStateless<DisplacementDim>& current_stateless,
              PrevState<SaturationData> const& S_L_prev_data,
              SaturationData& S_L_data, SaturationDataDeriv& dS_L_data) const;

private:
    SolidConstitutiveRelation<DisplacementDim> const& solid_material_;

    MSM::OGSMFrontTangentOperatorBlocksView<
        DisplacementDim,
        MSM::ForcesGradsCombinations<
            boost::mp11::mp_list<MSM::Strain, MSM::LiquidPressure>,
            boost::mp11::mp_list<MSM::Stress, MSM::Saturation>,
            boost::mp11::mp_list<MSM::Temperature>>::type>
        tangent_operator_blocks_view_;
};

extern template struct SolidMechanicsModel<2>;
extern template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
