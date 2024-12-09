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
std::pair<std::string, std::size_t> ParentDataType2String(ParentDataType p)
{
    // not used in OGS ParentDataType::POLYGON, ParentDataType::POLYHEDRON,
    // ParentDataType::HEXAHEDRON_24
    if (p == ParentDataType::MIXED)
    {
        return {"Mixed", 1};
    }

    if (p == ParentDataType::POLYVERTEX)
    {
        return {"Polyvertex", 1};
    }
    if (p == ParentDataType::POLYLINE)
    {
        return {"Polyline", 2};
    }
    if (p == ParentDataType::TRIANGLE)
    {
        return {"Triangle", 3};
    }
    if (p == ParentDataType::QUADRILATERAL)
    {
        return {"Quadrilateral", 4};
    }
    if (p == ParentDataType::TETRAHEDRON)
    {
        return {"Tetrahedron", 4};
    }
    if (p == ParentDataType::PYRAMID)
    {
        return {"Pyramid", 5};
    }
    if (p == ParentDataType::WEDGE)
    {
        return {"Wedge", 6};
    }
    if (p == ParentDataType::HEXAHEDRON)
    {
        return {"Hexahedron", 8};
    }
    if (p == ParentDataType::EDGE_3)
    {
        return {"Edge_3", 3};
    }
    if (p == ParentDataType::QUADRILATERAL_9)
    {
        return {"Quadrilateral_9", 9};
    }
    if (p == ParentDataType::TRIANGLE_6)
    {
        return {"Triangle_6", 6};
    }
    if (p == ParentDataType::QUADRILATERAL_8)
    {
        return {"Quadrilateral_8", 8};
    }
    if (p == ParentDataType::TETRAHEDRON_10)
    {
        return {"Tetrahedron_10", 10};
    }
    if (p == ParentDataType::PYRAMID_13)
    {
        return {"Pyramid_13", 13};
    }
    if (p == ParentDataType::WEDGE_15)
    {
        return {"Wedge_15", 15};
    }
    if (p == ParentDataType::WEDGE_18)
    {
        return {"Wedge_18", 18};
    }
    if (p == ParentDataType::HEXAHEDRON_20)
    {
        return {"Hexahedron_20", 20};
    }
    if (p == ParentDataType::HEXAHEDRON_27)
    {
        return {"Hexahedron_27", 27};
    }
    return {"Mixed", 1};
}
