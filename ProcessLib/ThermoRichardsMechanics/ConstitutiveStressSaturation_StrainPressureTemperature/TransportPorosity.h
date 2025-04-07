/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Bishops.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Porosity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SolidCompressibilityData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TransportPorosityData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{

template <int DisplacementDim>
struct TransportPorosityModel
{
    void eval(
        SpaceTimeData const& x_t, MediaData const& media_data,
        SolidCompressibilityData const& solid_compressibility_data,
        BishopsData const& bishops_data,
        PrevState<BishopsData> const& bishops_data_prev,
        CapillaryPressureData<DisplacementDim> const& p_cap_data,
        PorosityData const& poro_data,
        // TODO this should be mechanical strain, see the other TRM subtype
        StrainData<DisplacementDim> const& eps_data,
        PrevState<StrainData<DisplacementDim>> const& eps_prev_data,
        PrevState<TransportPorosityData> const& transport_poro_data_prev,
        TransportPorosityData& transport_poro_data) const;
};

extern template struct TransportPorosityModel<2>;
extern template struct TransportPorosityModel<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
