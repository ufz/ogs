// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TotalStress.h"

#include "MathLib/KelvinVector.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

template <int DisplacementDim>
void TotalStressModel<DisplacementDim>::eval(
    ProcessLib::ConstitutiveRelations::EffectiveStressData<
        DisplacementDim> const& eff_stress_data,
    BiotData const& biot_data,
    BishopsData const& chi_S_L,
    GasPressureData const& p_GR,
    CapillaryPressureData const& p_cap,
    TotalStressData<DisplacementDim>& out) const
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    out.sigma_total = (eff_stress_data.sigma_eff -
                       biot_data() * (p_GR.pG - chi_S_L.chi_S_L * p_cap.pCap) *
                           Invariants::identity2)
                          .eval();
}
template struct TotalStressModel<2>;
template struct TotalStressModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
