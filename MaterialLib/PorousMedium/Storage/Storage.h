/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   Storage.h
 *
 * Created on August 16, 2016, 12:53 PM
 */

#ifndef OGS_STORAGE_H
#define OGS_STORAGE_H

#include <string>

namespace MaterialLib
{
namespace PorousMedium
{
class Storage
{
public:

    virtual ~Storage()
    {
    }

    /// Get model name.
    virtual std::string getName() const = 0;

    /// Get property value.
    /// The argument is a double variable
    virtual double getValue(const double /*variable*/) const = 0;
};

}  // end of namespace
}  // end of namespace

#endif /* STORAGE_H */
