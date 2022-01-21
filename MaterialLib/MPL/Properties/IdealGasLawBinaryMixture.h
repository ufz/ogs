/**
 * \file
 * \author Norbert Grunwald
 * \date   Jan, 2021
 * \brief
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;
/**
 * \brief Density function for binary ideal gases.
 * \details This property must be a phase property, it
 * computes the density of a mixture of ideal gases as function of phase
 * pressure, capillary pressure, and temperature.
 */
class IdealGasLawBinaryMixture final : public Property
{
public:
    explicit IdealGasLawBinaryMixture(std::string name);

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'IdealGasLawBinaryMixture' is implemented on the "
                "'phase' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property values_ and dValues_.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/,
                             double const /*dt*/) const override;
};

}  // namespace MaterialPropertyLib
