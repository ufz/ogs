/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Advection.h"
#include "Base.h"
#include "Biot.h"
#include "ConstitutiveDensity.h"
#include "FluidDensity.h"
#include "PermeabilityData.h"
#include "PhaseTransitionData.h"
#include "Porosity.h"
#include "Saturation.h"
#include "SolidCompressibility.h"
#include "Viscosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct FW1Data
{
    GlobalDimMatrix<DisplacementDim> A;
};

template <int DisplacementDim>
struct FW1Model
{
    void eval(AdvectionData<DisplacementDim> const& advection_data,
              FluidDensityData const& fluid_density_data,
              FW1Data<DisplacementDim>& fW_1) const;
};

extern template struct FW1Model<2>;
extern template struct FW1Model<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
