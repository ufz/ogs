/**
 * \file
 * \author  Karsten Rink
 * \date    2012-05-02
 * \brief   Definition of mesh-related Enumerations.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

namespace MeshLib {

/**
 * \brief Types of mesh elements supported by OpenGeoSys.
 * Values are from VTKCellType enum
 */
enum class MeshElemType
{
    INVALID = 0,
    POINT = 1,
    LINE = 3,
    QUAD = 9,
    HEXAHEDRON = 12,
    TRIANGLE = 5,
    TETRAHEDRON = 10,
    PRISM = 16,
    PYRAMID = 14
};

/**
 * \brief Types of mesh elements supported by OpenGeoSys.
 */
enum class CellType
{
    INVALID = 0,
    POINT1 = 1,
    LINE2 = 2,
    LINE3 = 3,
    TRI3 = 4,
    TRI6 = 5,
    QUAD4 = 6,
    QUAD8 = 7,
    QUAD9 = 8,
    TET4 = 9,
    TET10 = 10,
    HEX8 = 11,
    HEX20 = 12,
    HEX27 = 13,
    PRISM6 = 14,
    PRISM15 = 15,
    PRISM18 = 16,
    PYRAMID5 = 17,
    PYRAMID13 = 18
};

/**
 * \brief Describes a mesh quality metric.
 */
enum class MeshQualityType
{
    INVALID = 0,
    ELEMENTSIZE,
    SIZEDIFFERENCE,
    EDGERATIO,
    EQUIANGLESKEW,
    RADIUSEDGERATIO
};

/**
 * \brief Selection of possible interpretations for intensities.
 */
enum class UseIntensityAs
{
    ELEVATION,
    MATERIALS,
    DATAVECTOR,
    NONE
};

/// Given a MeshElemType this returns the appropriate string.
const std::string MeshElemType2String(const MeshElemType t);

/// Given a MeshElemType this returns the appropriate string with a short name.
const std::string MeshElemType2StringShort(const MeshElemType t);

/// Given a string of the shortened name of the element type, this returns the corresponding MeshElemType.
MeshElemType String2MeshElemType(const std::string &s);

/// Returns a vector of all mesh element types
std::vector<MeshElemType> getMeshElemTypes();

/// Returns a vector of strings of mesh element types
std::vector<std::string> getMeshElemTypeStringsShort();

/// Given a MeshElemType this returns the appropriate string.
const std::string CellType2String(const CellType t);

const std::string MeshQualityType2String(const MeshQualityType t);

}
