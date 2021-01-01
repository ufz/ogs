/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
 * \brief Kozeny-Carman equation
 *
 * \f[
 * k = k_0 \left( \frac{1 - \phi_0}{1 - \phi} \right)^2 \left(
 * \frac{\phi}{\phi_0} \right)^3, \f] where \f$k\f$ is the permeability,
 * \fk_0\f$ is the initial permeability,
 * \f$\phi\f$ is the porosity, and
 * \f$\phi_0\f$ is the initial porosity.
 */
class KozenyCarmanModel final : public Property
{
public:
    explicit KozenyCarmanModel(ParameterLib::Parameter<double> const& k0,
                               ParameterLib::Parameter<double> const& phi0)
        : _k0(k0), _phi0(phi0)
    {
    }

    PropertyDataType value(
        MaterialPropertyLib::VariableArray const& variable_array,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const /*dt*/) const override;

private:
    /// Initial medium permeability
    ParameterLib::Parameter<double> const& _k0;
    /// Initial porosity
    ParameterLib::Parameter<double> const& _phi0;
};
}  // namespace MaterialPropertyLib
