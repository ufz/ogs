/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   Storage.h
 *
 * Created on August 16, 2016, 12:53 PM
 */

#pragma once

#include <string>

namespace MaterialLib
{
namespace PorousMedium
{
class Storage
{
public:
    virtual ~Storage() = default;
    /// Get model name.
    virtual std::string getName() const = 0;

    /**
     *  Get property value.
     *  @param variable A double type variable
     */
    virtual double getValue(const double variable) const = 0;
};

}  // end of namespace
}  // end of namespace
