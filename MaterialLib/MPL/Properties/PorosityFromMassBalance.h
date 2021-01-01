/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/// Porosity depending on the volumetric strain rate and phase pressure rate.
///
/// This property must be a solid phase.
class PorosityFromMassBalance final : public Property
{
private:
    /// Parameter, which is used by FEM to set the initial porosity value.
    ParameterLib::Parameter<double> const& phi0_;
    double const phi_min_;  //< Lower limit for the porosity.
    double const phi_max_;  //< Upper limit for the porosity.

public:
    PorosityFromMassBalance(
        std::string name,
        ParameterLib::Parameter<double> const& initial_porosity,
        double const phi_min, double const phi_max)
        : phi0_(initial_porosity), phi_min_(phi_min), phi_max_(phi_max)
    {
        name_ = std::move(name);
    }

    void checkScale() const override;

    PropertyDataType initialValue(ParameterLib::SpatialPosition const& pos,
                                  double const t) const override
    {
        return fromVector(phi0_(t, pos));
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType value(VariableArray const& variable_array,
                           VariableArray const& variable_array_prev,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
