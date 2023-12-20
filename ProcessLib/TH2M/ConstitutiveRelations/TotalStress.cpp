/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TotalStress.h"

#include "MathLib/KelvinVector.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

template <int DisplacementDim>
void TotalStressModel<DisplacementDim>::eval(
    ProcessLib::ConstitutiveRelations::StressData<DisplacementDim> const&
        eff_stress_data,
    BiotData const& biot_data,
    BishopsData const& chi_S_L,
    GasPressureData const& p_GR,
    CapillaryPressureData const& p_cap,
    TotalStressData<DisplacementDim>& out) const
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    out.sigma_total = (eff_stress_data.sigma -
                       biot_data() * (p_GR() - chi_S_L.chi_S_L * p_cap()) *
                           Invariants::identity2)
                          .eval();
}
template struct TotalStressModel<2>;
template struct TotalStressModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
