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

#include "Base.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct ElasticTangentStiffnessData
{
    KelvinMatrix<DisplacementDim> C_el;
};

template <int DisplacementDim>
struct ElasticTangentStiffnessModel
{
    explicit ElasticTangentStiffnessModel(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material_(solid_material)
    {
    }

    void eval(SpaceTimeData const& x_t,
              TemperatureData<DisplacementDim> const& T_data,
              ElasticTangentStiffnessData<DisplacementDim>& out) const;

private:
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;
};

extern template struct ElasticTangentStiffnessModel<2>;
extern template struct ElasticTangentStiffnessModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
