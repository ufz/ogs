/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Biot.h"
#include "ElasticTangentStiffnessData.h"
#include "Saturation.h"
#include "SolidCompressibility.h"
#include "SolidMechanics.h"
#include "SolidThermalExpansion.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
/// Data whose state must be tracked by the process.
template <int DisplacementDim>
struct StatefulData
{
    SaturationData S_L_data;

    static auto reflect()
    {
        using Self = StatefulData<DisplacementDim>;

        return Reflection::reflectWithoutName(&Self::S_L_data);
    }
};

template <int DisplacementDim>
struct StatefulDataPrev
{
    PrevState<SaturationData> S_L_data;

    StatefulDataPrev<DisplacementDim>& operator=(
        StatefulData<DisplacementDim> const& state)
    {
        S_L_data = state.S_L_data;

        return *this;
    }
};

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
struct ConstitutiveData
{
    SolidMechanicsDataStateless<DisplacementDim> s_mech_data;
};

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
struct ConstitutiveTempData
{
    ElasticTangentStiffnessData<DisplacementDim> C_el_data;
    BiotData biot_data;
    SolidCompressibilityData beta_p_SR;
    SaturationDataDeriv dS_L_dp_cap;
    SolidThermalExpansionData<DisplacementDim> s_therm_exp_data;

    using DisplacementDimVector = Eigen::Matrix<double, DisplacementDim, 1>;
    using DisplacementDimMatrix =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>;

    DisplacementDimMatrix dlambda_dp_GR;
    DisplacementDimMatrix dlambda_dp_cap;
    DisplacementDimMatrix dlambda_dT;
    DisplacementDimVector drho_GR_h_w_eff_dp_GR_Npart;
    DisplacementDimMatrix drho_GR_h_w_eff_dp_GR_gradNpart;
    DisplacementDimVector drho_LR_h_w_eff_dp_cap_Npart;
    DisplacementDimMatrix drho_LR_h_w_eff_dp_cap_gradNpart;
    DisplacementDimVector drho_GR_h_w_eff_dT;
    DisplacementDimMatrix dfW_4_LWpG_a_dp_GR;
    DisplacementDimMatrix dfW_4_LWpG_a_dp_cap;
    DisplacementDimMatrix dfW_4_LWpG_a_dT;
    DisplacementDimMatrix dfW_4_LWpG_d_dp_GR;
    DisplacementDimMatrix dfW_4_LWpG_d_dp_cap;
    DisplacementDimMatrix dfW_4_LWpG_d_dT;
    DisplacementDimMatrix dfW_4_LWpC_a_dp_GR;
    DisplacementDimMatrix dfW_4_LWpC_a_dp_cap;
    DisplacementDimMatrix dfW_4_LWpC_a_dT;
    DisplacementDimMatrix dfW_4_LWpC_d_dp_GR;
    DisplacementDimMatrix dfW_4_LWpC_d_dp_cap;
    DisplacementDimMatrix dfW_4_LWpC_d_dT;
    DisplacementDimMatrix dfC_4_LCpG_dT;
    DisplacementDimMatrix dfC_4_LCpC_a_dp_GR;
    DisplacementDimMatrix dfC_4_LCpC_a_dp_cap;
    DisplacementDimMatrix dfC_4_LCpC_a_dT;
    DisplacementDimMatrix dfC_4_LCpC_d_dp_GR;
    DisplacementDimMatrix dfC_4_LCpC_d_dp_cap;
    DisplacementDimMatrix dfC_4_LCpC_d_dT;
    DisplacementDimMatrix dadvection_C_dp_GR;
    DisplacementDimMatrix dadvection_C_dp_cap;
    DisplacementDimMatrix dk_over_mu_G_dp_cap;
    DisplacementDimMatrix dk_over_mu_L_dp_cap;
    double drho_u_eff_dT = std::numeric_limits<double>::quiet_NaN();
    double drho_u_eff_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double drho_u_eff_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double drho_h_eff_dT = std::numeric_limits<double>::quiet_NaN();
    double drho_h_eff_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double drho_h_eff_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dh_G_dT = std::numeric_limits<double>::quiet_NaN();
    double dh_L_dT = std::numeric_limits<double>::quiet_NaN();
    double dfC_4_MCpG_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfC_4_MCpG_dT = std::numeric_limits<double>::quiet_NaN();
    double dfC_4_MCT_dT = std::numeric_limits<double>::quiet_NaN();
    double dfC_4_MCu_dT = std::numeric_limits<double>::quiet_NaN();
    double dfC_3a_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfC_3a_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfC_3a_dT = std::numeric_limits<double>::quiet_NaN();
    double dfC_2a_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfC_2a_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfC_2a_dT = std::numeric_limits<double>::quiet_NaN();
    double dfW_2a_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfW_2b_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfW_2a_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfW_2b_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfW_2a_dT = std::numeric_limits<double>::quiet_NaN();
    double dfW_2b_dT = std::numeric_limits<double>::quiet_NaN();
    double dfW_3a_dp_GR = std::numeric_limits<double>::quiet_NaN();
    double dfW_3a_dp_cap = std::numeric_limits<double>::quiet_NaN();
    double dfW_3a_dT = std::numeric_limits<double>::quiet_NaN();
    double chi_s_L = std::numeric_limits<double>::quiet_NaN();
    double dchi_ds_L = std::numeric_limits<double>::quiet_NaN();
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
