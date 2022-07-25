/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Biot.h"
#include "ElasticTangentStiffness.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct SolidCompressibilityData
{
    double beta_SR;
};

template <int DisplacementDim>
struct SolidCompressibilityModel
{
    explicit SolidCompressibilityModel(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material_(solid_material)
    {
    }

    void eval(SpaceTimeData const& x_t,
              BiotData const& biot_data,
              ElasticTangentStiffnessData<DisplacementDim> const& C_el_data,
              SolidCompressibilityData& out) const;

private:
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;
};

extern template struct SolidCompressibilityModel<2>;
extern template struct SolidCompressibilityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
