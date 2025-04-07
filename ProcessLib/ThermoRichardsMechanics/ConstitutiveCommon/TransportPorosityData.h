/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib::ThermoRichardsMechanics
{

struct TransportPorosityData
{
    double phi;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName(
            "transport_porosity", &TransportPorosityData::phi);
    }
};

}  // namespace ProcessLib::ThermoRichardsMechanics
