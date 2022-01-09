/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
/**
 * Temperature dependent model for molecular diffusion in the pore water
 *
 * This model takes the form of
 * \f[
 * \mathbf{D} = \mathbf{D}_0 \mathrm{exp} \left( \frac{E_{\mathrm{a}}}{R}
 * \left( \frac{1}{T_0} - \frac{1}{T} \right)  \right),
 * \f]
 * where
 * \f$\mathbf{D_0}\f$ is the molecular diffusion at the reference temperature,
 * \f$E_{\mathrm{a}}\f$ is the activition energy for diffusion,
 * \f$R\f$ is the ideal gas constant,
 * \f$T\f$ is the absolute temperature,
 * \f$T_0\f$ is the reference temperature.
 * */
class TemperatureDependentDiffusion final : public Property
{
public:
    explicit TemperatureDependentDiffusion(
        ParameterLib::Parameter<double> const& D0,
        double const Ea,
        double const T0)
        : D0_(D0), Ea_(Ea), T0_(T0)
    {
    }

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const /*dt*/) const override;

private:
    /// the molecular diffusion at the reference temperature
    ParameterLib::Parameter<double> const& D0_;
    /// the activition energy for diffusion
    double const Ea_;
    /// the reference temperature
    double const T0_;
};
}  // namespace MaterialPropertyLib
