/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

template <int DisplacementDim>
struct ThermalConductivityData
{
    GlobalDimMatrix<DisplacementDim> lambda;
    // Currently unused, but there is a comment in TH2MFEM-impl.h referring to
    // this matrix
    // GlobalDimMatrix<DisplacementDim> dlambda_dp_GR;
    GlobalDimMatrix<DisplacementDim> dlambda_dp_cap;
    GlobalDimMatrix<DisplacementDim> dlambda_dT;
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
