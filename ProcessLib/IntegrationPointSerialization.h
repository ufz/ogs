/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{

/// Interface for serialization of process specific integration point data. The
/// interface requires a write and a read functions to be implemented in such a
/// way that so called round-trip is possible.
struct IntegrationPointSerialization
{
    virtual std::size_t writeIntegrationPointData(std::vector<char>& data) = 0;
    virtual void readIntegrationPointData(std::vector<char> const& data) = 0;

    virtual ~IntegrationPointSerialization() = default;
};

}   // namespace ProcessLib
