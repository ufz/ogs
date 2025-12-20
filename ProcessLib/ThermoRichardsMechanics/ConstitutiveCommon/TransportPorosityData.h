// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
