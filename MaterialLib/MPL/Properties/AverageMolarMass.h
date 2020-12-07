/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 07 2020
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * Molar mass of a mixture by mole fraction weighted average \f$ M=\Sigma_\zeta
 * x_{n,\alpha}^{\zeta}M^{\zeta}\f$.
 *
 * This property must be a phase property.
 */
class AverageMolarMass final : public Property
{
public:
    explicit AverageMolarMass(std::string name);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
