/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MathLib/KelvinVector.h"

namespace ProcessLib::TH2M
{
/// Variables needed only for the assembly process. The values are not preserved
/// throughout the iterations contrary to the variables in IntegrationPointData.
template <int DisplacementDim>
struct ConstitutiveVariables
{
    using KelvinMatrixType =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using DisplacementDimVector = Eigen::Matrix<double, DisplacementDim, 1>;
    using DisplacementDimMatrix =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>;

    KelvinMatrixType C;

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
    DisplacementDimMatrix dfC_4_LCpG_dT;
    DisplacementDimMatrix dadvection_C_dp_GR;
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
    double ds_L_dp_cap = std::numeric_limits<double>::quiet_NaN();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib::TH2M
