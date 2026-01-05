// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

namespace FileIO
{
namespace Gocad
{
enum class DataType
{
    UNDEFINED,
    VSET,
    PLINE,
    TSURF,
    MODEL3D,
    ALL
};

/// Given a Gocad DataType this returns the appropriate string.
std::string dataType2String(DataType const t);

/// Given a Gocad DataType this returns the appropriate short form.
std::string dataType2ShortString(DataType const t);

}  // namespace Gocad

}  // namespace FileIO
