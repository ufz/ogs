/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 16, 2016, 12:53 PM
 */

#pragma once

#include <string>

#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace PorousMedium
{
class Porosity
{
public:
    explicit Porosity(ParameterLib::Parameter<double> const& parameter)
        : _parameter(parameter)
    {
    }
    virtual ~Porosity() = default;

    /**
     *  Get property value.
     *  @param t point in time
     *  @param pos spatial position
     *  @param variable    A variable with any double type value.
     *  @param temperature Temperature with any double type value.
     */
    virtual double getValue(const double t,
                            ParameterLib::SpatialPosition const& pos,
                            const double variable,
                            const double temperature) const
    {
        (void)variable;
        (void)temperature;
        return _parameter(t, pos)[0];
    }

private:
    ParameterLib::Parameter<double> const& _parameter;
};

}  // namespace PorousMedium
}  // namespace MaterialLib
