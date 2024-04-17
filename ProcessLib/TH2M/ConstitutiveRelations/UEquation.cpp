/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "UEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
void FU2KUpCModel::eval(BiotData const& biot_data,
                        BishopsData const& chi_S_L,
                        FU2KUpCData& fu_2_KupC) const
{
    fu_2_KupC.m = biot_data() * chi_S_L.chi_S_L;
}

void FU2KUpCModel::dEval(BiotData const& biot_data,
                         BishopsData const& chi_S_L,
                         CapillaryPressureData const& p_cap,
                         SaturationDataDeriv const& dS_L_dp_cap,
                         FU2KUpCDerivativeData& dfu_2_KupC) const
{
    dfu_2_KupC.dp_cap =
        biot_data() * chi_S_L.dchi_dS_L * dS_L_dp_cap() * p_cap();
}

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
