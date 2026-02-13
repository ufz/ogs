// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace ProcessLib::ThermoRichardsMechanics
{
struct TRMStorageData
{
    double storage_p_a_p;
    double storage_p_a_S_X_NTN;
    double J_pp_X_NTN;
    double storage_p_a_S_Jpp_X_NTN;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
