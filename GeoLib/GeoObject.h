// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "GeoType.h"

namespace GeoLib
{
struct GeoObject
{
    virtual ~GeoObject() = default;
    /// return a geometry type
    virtual GEOTYPE getGeoType() const = 0;
};
}  // end namespace GeoLib
