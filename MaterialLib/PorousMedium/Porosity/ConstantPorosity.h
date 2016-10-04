/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   ConstantPorosity.h
 *
 * Created on August 16, 2016, 1:03 PM
 */

#ifndef OGS_CONSTANTPOROSITY_H
#define OGS_CONSTANTPOROSITY_H

#include "Porosity.h"

namespace MaterialLib
{
namespace PorousMedium
{
class ConstantPorosity final : public Porosity
{
public:

    explicit ConstantPorosity(const double value) : _value(value)
    {
    }

    /// Get model name.
    std::string getName() const override
    {
        return "Constant porosity";
    }

    /**
     *  Get property value.
     *  \param variable    A variable with any double type value.
     *  \param temperature Temperature with any double type value.
     */
    double getValue(const double variable, double temperature) const override
    {
        (void) variable;
        (void) temperature;
        return _value;
    }

private:
    const double _value;
};

}  // end of namespace
}  // end of namespace

#endif /* CONSTANTPOROSITY_H */
