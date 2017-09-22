/**
 * \author Norbert Grunwald
 * \date   11.09.2017
 * \brief  
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pLinearTemperature.h"
#include "../mpComponent.h"
#include <iostream>

namespace MaterialPropertyLib
{
/// This constructor throws an error since it was used for a medium
/// property while it is a component property.
LinearTemperature::LinearTemperature(Medium*)
: _component (0)
{
    notImplemented("LinearTemperature", "Medium");
};

/// This constructor throws an error since it was used for a phase
/// property while it is a component property.
LinearTemperature::LinearTemperature(Phase*)
: _component (0)
{
    notImplemented("LinearTemperature", "Phase");
};

/// This empty constructor is probably unnecessary.
LinearTemperature::LinearTemperature()
: _component (0)
{};
/// This constructor copies the pointer of the component from the
/// arguments into the private attribute.
LinearTemperature::LinearTemperature(Component* c)
: _component (c)
{};

/**
 * This method computes a density by a linear correlation with temperature
 * by \f$\rho=\rho_0+\partial\rho/\partial T\left(T-T_\mathrm{ref}\right)\f$
 */
PropertyDataType LinearTemperature::value(VariableArray const& v)
{
	const double rho_0 = getScalar(_component->property(reference_density));
	const double drho_dT  = getScalar(_component->property(drhodT));
	const double T_ref = getScalar(_component->property(reference_temperature));

	_value = rho_0 + drho_dT*(getScalar(v[T])-T_ref);
	return _value;
}

} //MaterialPropertyLib

