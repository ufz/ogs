// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Medium.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"

namespace ProcessLib::SmallDeformation
{

using namespace ProcessLib::ConstitutiveRelations;

struct MediaData
{
    explicit MediaData(MaterialPropertyLib::Medium const& medium)
        : medium{medium}, solid{medium.phase("Solid")}
    {
    }

    MaterialPropertyLib::Medium const& medium;
    MaterialPropertyLib::Phase const& solid;
};

}  // namespace ProcessLib::SmallDeformation
