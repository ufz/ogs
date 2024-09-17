/**
 * \file
 * \brief  Enum ParentDataType to string translations
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshPropertyDataType.h"

// See https://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology (Arbitrary)
std::string ParentDataType2String(ParentDataType p)
{
    // not used in OGS ParentDataType::POLYGON, ParentDataType::POLYHEDRON,
    // ParentDataType::HEXAHEDRON_24
    if (p == ParentDataType::MIXED)
    {
        return "Mixed";
    }

    if (p == ParentDataType::POLYVERTEX)
    {
        return "Polyvertex";
    }
    if (p == ParentDataType::POLYLINE)
    {
        return "Polyline";
    }
    if (p == ParentDataType::TRIANGLE)
    {
        return "Triangle";
    }
    if (p == ParentDataType::QUADRILATERAL)
    {
        return "Quadrilateral";
    }
    if (p == ParentDataType::TETRAHEDRON)
    {
        return "Tetrahedron";
    }
    if (p == ParentDataType::PYRAMID)
    {
        return "Pyramid";
    }
    if (p == ParentDataType::WEDGE)
    {
        return "Wedge";
    }
    if (p == ParentDataType::HEXAHEDRON)
    {
        return "Hexahedron";
    }
    if (p == ParentDataType::EDGE_3)
    {
        return "Edge_3";
    }
    if (p == ParentDataType::QUADRILATERAL_9)
    {
        return "Quadrilateral_9";
    }
    if (p == ParentDataType::TRIANGLE_6)
    {
        return "Triangle_6";
    }
    if (p == ParentDataType::QUADRILATERAL_8)
    {
        return "Quadrilateral_8";
    }
    if (p == ParentDataType::TETRAHEDRON_10)
    {
        return "Tetrahedron_10";
    }
    if (p == ParentDataType::PYRAMID_13)
    {
        return "Pyramid_13";
    }
    if (p == ParentDataType::WEDGE_15)
    {
        return "Wedge_15";
    }
    if (p == ParentDataType::WEDGE_18)
    {
        return "Wedge_18";
    }
    if (p == ParentDataType::HEXAHEDRON_20)
    {
        return "Hexahedron_20";
    }
    if (p == ParentDataType::HEXAHEDRON_27)
    {
        return "Hexahedron_27";
    }
    return "Mixed";
}
