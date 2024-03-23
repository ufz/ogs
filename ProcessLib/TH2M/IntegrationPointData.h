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

#include <memory>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace TH2M
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;
    using GlobalDimVectorType =
        typename ShapeMatricesTypePressure::GlobalDimVectorType;

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    // phase intrinsic densities
    double rhoGR = std::numeric_limits<double>::quiet_NaN();
    double rhoLR = std::numeric_limits<double>::quiet_NaN();
    double rhoSR = std::numeric_limits<double>::quiet_NaN();

    // phase intrinsic density derivatives
    double drhoGR_dpGR = std::numeric_limits<double>::quiet_NaN();
    double drhoGR_dpCap = std::numeric_limits<double>::quiet_NaN();
    double drhoGR_dT = std::numeric_limits<double>::quiet_NaN();
    double drhoLR_dpGR = std::numeric_limits<double>::quiet_NaN();
    double drhoLR_dpCap = std::numeric_limits<double>::quiet_NaN();
    double drhoLR_dT = std::numeric_limits<double>::quiet_NaN();

    // vapour pressure (water component partial gas phase pressure)
    double pWGR = std::numeric_limits<double>::quiet_NaN();

    // real constituent partial densities
    double rhoCGR = std::numeric_limits<double>::quiet_NaN();
    double rhoCGR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoWGR = std::numeric_limits<double>::quiet_NaN();
    double rhoWGR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoCLR = std::numeric_limits<double>::quiet_NaN();
    double rhoCLR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoWLR = std::numeric_limits<double>::quiet_NaN();
    double rhoWLR_prev = std::numeric_limits<double>::quiet_NaN();

    // real constituent partial density derivatives
    double drhoCGR_dpGR = std::numeric_limits<double>::quiet_NaN();
    double drhoWGR_dpGR = std::numeric_limits<double>::quiet_NaN();
    double drhoCGR_dpCap = std::numeric_limits<double>::quiet_NaN();
    double drhoWGR_dpCap = std::numeric_limits<double>::quiet_NaN();
    double drhoCGR_dT = std::numeric_limits<double>::quiet_NaN();
    double drhoWGR_dT = std::numeric_limits<double>::quiet_NaN();

    // phase composition
    // molar fraction
    double xnCG = std::numeric_limits<double>::quiet_NaN();
    double xnWG = std::numeric_limits<double>::quiet_NaN();
    double xnWL = std::numeric_limits<double>::quiet_NaN();

    // mass fraction
    double xmCG = std::numeric_limits<double>::quiet_NaN();
    double xmWG = std::numeric_limits<double>::quiet_NaN();
    double xmWL = std::numeric_limits<double>::quiet_NaN();

    // mass fraction derivatives
    double dxmWG_dpGR = std::numeric_limits<double>::quiet_NaN();
    double dxmWG_dpCap = std::numeric_limits<double>::quiet_NaN();
    double dxmWG_dT = std::numeric_limits<double>::quiet_NaN();
    double dxmWL_dpGR = std::numeric_limits<double>::quiet_NaN();
    double dxmWL_dpCap = std::numeric_limits<double>::quiet_NaN();
    double dxmWL_dpLR = std::numeric_limits<double>::quiet_NaN();
    double dxmWL_dT = std::numeric_limits<double>::quiet_NaN();

    // diffusion coefficients
    double diffusion_coefficient_vapour =
        std::numeric_limits<double>::quiet_NaN();
    double diffusion_coefficient_solute =
        std::numeric_limits<double>::quiet_NaN();

    // phase enthalpies
    double rho_G_h_G = std::numeric_limits<double>::quiet_NaN();
    double rho_G_h_G_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S_prev = std::numeric_limits<double>::quiet_NaN();

    // specific enthalpies
    double h_G = std::numeric_limits<double>::quiet_NaN();
    double h_L = std::numeric_limits<double>::quiet_NaN();
    double h_S = std::numeric_limits<double>::quiet_NaN();
    double h_CG = std::numeric_limits<double>::quiet_NaN();
    double h_WG = std::numeric_limits<double>::quiet_NaN();
    double h_WL = std::numeric_limits<double>::quiet_NaN();
    double h_CL = std::numeric_limits<double>::quiet_NaN();

    // internal energies
    double rho_u_eff = std::numeric_limits<double>::quiet_NaN();
    double rho_u_eff_prev = std::numeric_limits<double>::quiet_NaN();

    // porosity
    double phi = std::numeric_limits<double>::quiet_NaN();
    double dphi_dT = std::numeric_limits<double>::quiet_NaN();

    double muGR = std::numeric_limits<double>::quiet_NaN();
    double muLR = std::numeric_limits<double>::quiet_NaN();

    GlobalDimMatrixType lambda;
    GlobalDimVectorType d_CG;
    GlobalDimVectorType d_WG;
    GlobalDimVectorType d_CL;
    GlobalDimVectorType d_WL;

    GlobalDimVectorType w_GS;
    GlobalDimVectorType w_LS;

    double thermal_volume_strain = std::numeric_limits<double>::quiet_NaN();

    double integration_weight = std::numeric_limits<double>::quiet_NaN();

    void pushBackState()
    {
        rho_G_h_G_prev = rho_G_h_G;
        rho_L_h_L_prev = rho_L_h_L;
        rho_S_h_S_prev = rho_S_h_S;

        rhoCGR_prev = rhoCGR;
        rhoWGR_prev = rhoWGR;
        rhoCLR_prev = rhoCLR;
        rhoWLR_prev = rhoWLR;

        rho_u_eff_prev = rho_u_eff;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH2M
}  // namespace ProcessLib
