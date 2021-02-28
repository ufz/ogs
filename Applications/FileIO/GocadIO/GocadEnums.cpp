/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GocadEnums.h"

namespace FileIO
{
namespace Gocad
{
std::string dataType2String(DataType const t)
{
    if (t == DataType::UNDEFINED)
    {
        return "[Undefined Type]";
    }
    if (t == DataType::VSET)
    {
        return "GOCAD VSet";
    }
    if (t == DataType::PLINE)
    {
        return "GOCAD PLine";
    }
    if (t == DataType::TSURF)
    {
        return "GOCAD TSurf";
    }
    if (t == DataType::MODEL3D)
    {
        return "GOCAD Model3d";
    }
    return "[all types]";
}

std::string dataType2ShortString(DataType const t)
{
    if (t == DataType::UNDEFINED)
    {
        return "[undefined]";
    }
    if (t == DataType::VSET)
    {
        return "vertices";
    }
    if (t == DataType::PLINE)
    {
        return "line";
    }
    if (t == DataType::TSURF)
    {
        return "surface";
    }
    if (t == DataType::MODEL3D)
    {
        return "model";
    }
    return "[all data]";
}

}  // namespace Gocad

}  // namespace FileIO
