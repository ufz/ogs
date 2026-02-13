// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Medium.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct MediaData
{
    explicit MediaData(MaterialPropertyLib::Medium const& medium)
        : medium{medium},
          liquid{medium.phase("AqueousLiquid")},
          solid{medium.phase("Solid")}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& liquid;
    MaterialPropertyLib::Phase const& solid;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
