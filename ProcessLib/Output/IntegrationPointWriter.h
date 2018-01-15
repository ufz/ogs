/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
struct IntegrationPointWriter
{
    virtual ~IntegrationPointWriter() = default;

    virtual int numberOfComponents() const = 0;
    virtual int integrationOrder() const = 0;
    virtual std::string name() const = 0;
    virtual std::vector<std::vector<double>> values() const = 0;
};

}  // namespace ProcessLib
