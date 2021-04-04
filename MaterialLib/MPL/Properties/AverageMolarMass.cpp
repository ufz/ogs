/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 07 2020
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/AverageMolarMass.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
AverageMolarMass::AverageMolarMass(std::string name)
{
    name_ = std::move(name);
}

void AverageMolarMass::checkScale() const
{
    if (!(std::holds_alternative<Phase*>(scale_)))
    {
        OGS_FATAL(
            "The property 'AverageMolarMass' is implemented on the 'phase' "
            "scale only.");
    }
}

PropertyDataType AverageMolarMass::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto phase = std::get<Phase*>(scale_);
    auto const numberOfComponents = phase->numberOfComponents();
    if (numberOfComponents < 1)
    {
        // if phase contains no components, jst return phase molar mass
        return phase->property(PropertyType::molar_mass)
            .template value<double>(variable_array, pos, t, dt);
    }
    else if (numberOfComponents > 2)
    {
        OGS_FATAL(
            "AverageMolarMass::value only allows for phases consisting of up "
            "to two components.");
    }

    // TODO (grunwald) : Task here is to retrieve the individual molar fractions
    // of each compontne in the phase. Those are not static properties (as in
    // case of molar mass), but they depend on some state-dependent rule. Option
    // 1 is to call that rule here, option 2 would be to wrap this composition
    // info into some container and to put it into the variable_array. This
    // would have to be a vector of variable size, depending on the number of
    // components. Unfortunately, the data types of MPL::VariableArray do not
    // include such a type.
    //
    // Therefore, I go with option 1 here, which unfortunately only allows the
    // use of binary mixtures at the moment.
    auto const molar_fraction =
        phase->property(PropertyType::mole_fraction)
            .template value<Eigen::Vector2d>(variable_array, pos, t, dt);

    double M = 0.;
    for (size_t c = 0; c < numberOfComponents; c++)
    {
        auto const M_zeta =
            phase->component(c)
                .property(PropertyType::molar_mass)
                .template value<double>(variable_array, pos, t, dt);
        auto const xn_zeta = molar_fraction[c];

        M += xn_zeta * M_zeta;
    }

    return M;
}

PropertyDataType AverageMolarMass::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)primary_variable;
    assert(((primary_variable == Variable::phase_pressure) ||
            (primary_variable == Variable::temperature)) &&
           "AverageMolarMass::dValue is implemented for derivatives with "
           "respect to phase_pressure or temperature only.");

    auto phase = std::get<Phase*>(scale_);

    auto const numberOfComponents = phase->numberOfComponents();
    if (numberOfComponents <= 1)
    {
        return 0.;
    }
    else if (numberOfComponents > 2)
    {
        OGS_FATAL(
            "AverageMolarMass::dvalue is currently implemented two or less "
            "phase components only.");
    }

    // TODO (grunwald) : This should return a vector of length
    // phase->numberOfComponents(). Currently, this feature is implemented
    // for binary phases only.
    auto const dxnC = phase->property(PropertyType::mole_fraction)
                          .template dValue<Eigen::Vector2d>(
                              variable_array, primary_variable, pos, t, dt)[0];

    auto const M_0 = phase->component(0)
                         .property(PropertyType::molar_mass)
                         .template value<double>(variable_array, pos, t, dt);
    auto const M_1 = phase->component(1)
                         .property(PropertyType::molar_mass)
                         .template value<double>(variable_array, pos, t, dt);

    return dxnC * (M_0 - M_1);

}  // namespace MaterialPropertyLib

PropertyDataType AverageMolarMass::d2Value(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable1*/, Variable const /*primary_variable2*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("AverageMolarMass::d2Value is not yet implemented.");

    return 0.;
}

}  // namespace MaterialPropertyLib
