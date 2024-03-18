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
    double rhoSR = std::numeric_limits<double>::quiet_NaN();

    // real constituent partial densities
    double rhoCGR = std::numeric_limits<double>::quiet_NaN();
    double rhoCGR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoWGR = std::numeric_limits<double>::quiet_NaN();
    double rhoWGR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoCLR = std::numeric_limits<double>::quiet_NaN();
    double rhoCLR_prev = std::numeric_limits<double>::quiet_NaN();
    double rhoWLR = std::numeric_limits<double>::quiet_NaN();
    double rhoWLR_prev = std::numeric_limits<double>::quiet_NaN();

    // phase enthalpies
    double rho_G_h_G = std::numeric_limits<double>::quiet_NaN();
    double rho_G_h_G_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L = std::numeric_limits<double>::quiet_NaN();
    double rho_L_h_L_prev = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S = std::numeric_limits<double>::quiet_NaN();
    double rho_S_h_S_prev = std::numeric_limits<double>::quiet_NaN();

    // specific enthalpies
    double h_S = std::numeric_limits<double>::quiet_NaN();

    // internal energies
    double rho_u_eff = std::numeric_limits<double>::quiet_NaN();
    double rho_u_eff_prev = std::numeric_limits<double>::quiet_NaN();

    // porosity
    double phi = std::numeric_limits<double>::quiet_NaN();
    double dphi_dT = std::numeric_limits<double>::quiet_NaN();

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
