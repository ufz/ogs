// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

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
    char_native,
    uchar,
    enum_length
};

// See https://www.xdmf.org/index.php/XDMF_Model_and_Format - only relevant if
// TopologyType(s) added, structure is open to be extended with GridType(s)
// and ItemType(s).
enum class ParentDataType
{
    MIXED = 0,
    POLYVERTEX = 1,
    POLYLINE = 2,  // OGS polylines are supposed to contain exactly 2 nodes
    // POLYGON = 3, // not used in OGS
    TRIANGLE = 4,
    QUADRILATERAL = 5,
    TETRAHEDRON = 6,
    PYRAMID = 7,
    WEDGE = 8,
    HEXAHEDRON = 9,
    // POLYHEDRON = 16,  // not used in OGS
    EDGE_3 = 34,
    QUADRILATERAL_9 = 35,
    TRIANGLE_6 = 36,
    QUADRILATERAL_8 = 37,
    TETRAHEDRON_10 = 38,
    PYRAMID_13 = 39,
    WEDGE_15 = 40,
    WEDGE_18 = 41,
    HEXAHEDRON_20 = 48,
    // HEXAHEDRON_24 = 49,  // not used in OGS
    HEXAHEDRON_27 = 50
};

std::pair<std::string, std::size_t> ParentDataType2String(ParentDataType p);
