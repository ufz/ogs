/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
std::string dataType2String(DataType const t);

/// Given a Gocad DataType this returns the appropriate short form.
std::string dataType2ShortString(DataType const t);

}  // namespace Gocad

}  // namespace FileIO
