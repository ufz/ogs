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

#ifndef CONSTANTPOROSITY_H
#define CONSTANTPOROSITY_H

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

    /// Get porosity value.
    /// \param vars Variable values
    double getValue(const double vars[]) const override
    {
        (void) vars;
        return _value;
    }

private:
    const double _value;
};

}  // end of namespace
}  // end of namespace

#endif /* CONSTANTPOROSITY_H */
