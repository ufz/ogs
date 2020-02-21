/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    Phase* _phase = nullptr;

    /// Parameter, which is used by FEM to set the initial porosity value.
    ParameterLib::Parameter<double> const& _phi0;

public:
    PorosityFromMassBalance(
        ParameterLib::Parameter<double> const& initial_porosity)
        : _phi0(initial_porosity)
    {
    }

    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override;

    PropertyDataType initialValue(ParameterLib::SpatialPosition const& pos,
                                  double const t) const override
    {
        return fromVector(_phi0(t, pos));
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
