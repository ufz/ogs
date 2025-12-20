// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/ConstitutiveRelations/MechanicalStrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Bishops.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Porosity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidCompressibilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TransportPorosityData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{

template <int DisplacementDim>
struct TransportPorosityModel
{
    void eval(SpaceTimeData const& x_t, MediaData const& media_data,
              SolidCompressibilityData const& solid_compressibility_data,
              BishopsData const& bishops_data,
              PrevState<BishopsData> const& bishops_data_prev,
              CapillaryPressureData<DisplacementDim> const& p_cap_data,
              PorosityData const& poro_data,
              ProcessLib::ConstitutiveRelations::MechanicalStrainData<
                  DisplacementDim> const& eps_m_data,
              PrevState<ProcessLib::ConstitutiveRelations::MechanicalStrainData<
                  DisplacementDim>> const& eps_m_prev_data,
              PrevState<TransportPorosityData> const& transport_poro_data_prev,
              TransportPorosityData& transport_poro_data) const;
};

extern template struct TransportPorosityModel<2>;
extern template struct TransportPorosityModel<3>;
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
