// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"
#include "Biot.h"
#include "Bishops.h"
#include "ProcessLib/ConstitutiveRelations/EffectiveStressData.h"
#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

template <int DisplacementDim>
struct TotalStressData
{
    KelvinVector<DisplacementDim> sigma_total = KV::KVnan<DisplacementDim>();

    static auto reflect()
    {
        using Self = TotalStressData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("sigma_total",
                                                       &Self::sigma_total);
    }
};

template <int DisplacementDim>
struct TotalStressModel
{
    void eval(ProcessLib::ConstitutiveRelations::EffectiveStressData<
                  DisplacementDim> const& eff_stress_data,
              BiotData const& biot_data,
              BishopsData const& chi_S_L,
              GasPressureData const& p_GR,
              CapillaryPressureData const& p_cap,
              TotalStressData<DisplacementDim>& out) const;
};

extern template struct TotalStressModel<2>;
extern template struct TotalStressModel<3>;

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
