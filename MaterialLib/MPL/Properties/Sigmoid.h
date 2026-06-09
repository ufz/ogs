// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/**
 * Sigmoid property with configurable independent variable
 *
 * \details The sigmoid function is a smooth step function
 * with a configurable range defined by the following formula
 *
 * \f[
 *      (X_\mathrm{u} - X_\mathrm{l}) \left[1 + \exp(-k(X - X_\mathrm{c}))
 * \right]^{-1} + X_\mathrm{l}
 * \f]
 *
 * where
 * - \f$X_\mathrm{l}\f$ is the lower bound (value as X -> -inf)
 * - \f$X_\mathrm{u}\f$ is the upper bound (value as X -> +inf)
 * - \f$X\f$ is the independent material variable (temperature, pressure, etc.)
 * - \f$X_\mathrm{c}\f$ is the midpoint (a.k.a. inflection point) of the
 *   transition
 * - \f$k\f$ controls the steepness of the transition; a positive \f$k\f$ gives
 *   a monotonically increasing curve, a negative \f$k\f$ a decreasing one.
 *
 */
class Sigmoid final : public Property
{
public:
    Sigmoid(std::string name,
            double const steepness,
            double const characteristic_value,
            double const lower_bound,
            double const upper_bound,
            Variable const independent_variable);

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1,
                             Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t,
                             double const dt) const override;

private:
    /// True when the derivative w.r.t. \c variable is identically zero, i.e.
    /// \c variable is not the independent variable.
    bool derivativeVanishes(Variable const variable) const;

    double steepness_;  ///< \f$k\f$ parameter (controls steepness); positive k
                        ///< gives an increasing curve, negative k a decreasing
                        ///< one. Has units of \f$1/[X]\f$, i.e. the inverse of
                        ///< the independent variable's unit.
    double characteristic_value_;    ///< X_c parameter, the midpoint (a.k.a.
                                     ///< inflection point) of the transition
    double lower_bound_;             ///< X_l parameter (lower bound)
    double upper_bound_;             ///< X_u parameter (upper bound)
    Variable independent_variable_;  ///< Independent material variable
                                     ///< (temperature, pressure, etc.)
};

}  // namespace MaterialPropertyLib
