/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   Porosity.h
 *
 * Created on August 16, 2016, 12:53 PM
 */

#ifndef POROSITY_H
#define POROSITY_H

#include <string>

namespace MaterialLib
{
namespace PorousMedium
{
class Porosity
{
public:

    virtual ~Porosity()
    {
    }

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get property value.
    /// The argument is an array of variables.
    virtual double getValue(const double /* vars*/[] = nullptr) const = 0;
};

}  // end of namespace
}  // end of namespace

#endif /* POROSITY_H */
