/**
 * \file
 * \author Tobias Meisel
 * \date   2020-12-15
 * \brief  Enum for all propertyVector data types
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

// TODO (tm) If used on several other places move definition of propertyVector
enum class MeshPropertyDataType
{
    unknown = 0,
    float64,
    float32,
    int32,
    int64,
    uint32,
    uint64,
    int8,
    uint8,
    enum_length
};