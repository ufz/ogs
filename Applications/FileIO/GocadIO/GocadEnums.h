/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
std::string DataType2Str(DataType const t);

/// Given a Gocad DataType this returns the appropriate short form.
std::string DataType2ShortStr(DataType const t);

}  // namespace GocadIO

}  // namespace FileIO
