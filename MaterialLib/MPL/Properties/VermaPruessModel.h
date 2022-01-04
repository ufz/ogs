/**
 * \file
 *
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
 * \brief Verma-Pruess equation \cite verma1988thermohydrological
 *
 * \f[
 * k = k_0 \left( \frac{\phi - \phi_c}{\phi_0 - \phi_c} \right)^n, \f]
 * where \f$k\f$ is the permeability,
 * \f$k_0\f$ is the initial permeability,
 * \f$\phi\f$ is the porosity,
 * \f$\phi_0\f$ is the initial porosity,
 * \f$\phi_c\f$ is the critical porosity, and
 * \f$n\f$ is the exponent.
 */
class VermaPruessModel final : public Property
{
public:
    explicit VermaPruessModel(ParameterLib::Parameter<double> const& k0,
                              ParameterLib::Parameter<double> const& phi0,
                              ParameterLib::Parameter<double> const& phi_c,
                              ParameterLib::Parameter<double> const& n)
        : _k0(k0), _phi0(phi0), _phi_c(phi_c), _n(n)
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
    /// Critical porosity
    ParameterLib::Parameter<double> const& _phi_c;
    /// Exponent
    ParameterLib::Parameter<double> const& _n;
};
}  // namespace MaterialPropertyLib
