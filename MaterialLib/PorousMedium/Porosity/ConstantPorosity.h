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
class ConstantPorosity : public Porosity
{
public:
    ConstantPorosity(const double value) : _value(value) {}
    /// Get model name.
    virtual std::string getName() const final { return "Constant porosity"; }
    /// Get porosity value.
    /// \param vars Variable values
    virtual double getValue(const double /* vars*/[] = nullptr) const final
    {
        return _value;
    }

private:
    double _value;
};

}  // end of namespace
}  // end of namespace

#endif /* CONSTANTPOROSITY_H */
