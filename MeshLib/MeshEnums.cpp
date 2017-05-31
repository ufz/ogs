/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of mesh-related enumerations.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshEnums.h"

namespace MeshLib {

const std::string MeshElemType2String(const MeshElemType t)
{
    if (t == MeshElemType::POINT)
        return "Point";
    if (t == MeshElemType::LINE)
        return "Line";
    if (t == MeshElemType::QUAD)
        return "Quad";
    if (t == MeshElemType::HEXAHEDRON)
        return "Hexahedron";
    if (t == MeshElemType::TRIANGLE)
        return "Triangle";
    if (t == MeshElemType::TETRAHEDRON)
        return "Tetrahedron";
    if (t == MeshElemType::PRISM)
        return "Prism";
    if (t == MeshElemType::PYRAMID)
        return "Pyramid";
    return "none";
}

const std::string MeshElemType2StringShort(const MeshElemType t)
{
    if (t == MeshElemType::POINT)
        return "point";
    if (t == MeshElemType::LINE)
        return "line";
    if (t == MeshElemType::QUAD)
        return "quad";
    if (t == MeshElemType::HEXAHEDRON)
        return "hex";
    if (t == MeshElemType::TRIANGLE)
        return "tri";
    if (t == MeshElemType::TETRAHEDRON)
        return "tet";
    if (t == MeshElemType::PRISM)
        return "pris";
    if (t == MeshElemType::PYRAMID)
        return "pyra";
    return "none";
}

MeshElemType String2MeshElemType(const std::string &s)
{
    if ((s == "point") || (s == "Point"))
        return MeshElemType::POINT;
    if ((s == "line") || (s == "Line"))
        return MeshElemType::LINE;
    if ((s == "quad") || (s == "Quadrilateral"))
        return MeshElemType::QUAD;
    if ((s == "hex") || (s == "Hexahedron"))
        return MeshElemType::HEXAHEDRON;
    if ((s == "tri") || (s == "Triangle"))
        return MeshElemType::TRIANGLE;
    if ((s == "tet") || (s == "Tetrahedron"))
        return MeshElemType::TETRAHEDRON;
    if ((s == "pris") || (s == "Prism"))
        return MeshElemType::PRISM;
    if ((s == "pyra") || (s == "Pyramid"))
        return MeshElemType::PYRAMID;
    return MeshElemType::INVALID;
}

std::vector<MeshElemType> getMeshElemTypes()
{
    std::vector<MeshElemType> vec;
    vec.push_back(MeshElemType::POINT);
    vec.push_back(MeshElemType::LINE);
    vec.push_back(MeshElemType::QUAD);
    vec.push_back(MeshElemType::HEXAHEDRON);
    vec.push_back(MeshElemType::TRIANGLE);
    vec.push_back(MeshElemType::TETRAHEDRON);
    vec.push_back(MeshElemType::PRISM);
    vec.push_back(MeshElemType::PYRAMID);
    return vec;
}

std::vector<std::string> getMeshElemTypeStringsShort()
{
    std::vector<std::string> vec;
    for (MeshElemType eleType : getMeshElemTypes())
        vec.push_back(MeshElemType2StringShort(eleType));
    return vec;
}

const std::string CellType2String(const CellType t)
{
#define RETURN_CELL_TYPE_STR(t, type)\
    if ((t) == CellType::type)\
        return #type;

    RETURN_CELL_TYPE_STR(t, POINT1);
    RETURN_CELL_TYPE_STR(t, LINE2);
    RETURN_CELL_TYPE_STR(t, LINE3);
    RETURN_CELL_TYPE_STR(t, QUAD4);
    RETURN_CELL_TYPE_STR(t, QUAD8);
    RETURN_CELL_TYPE_STR(t, QUAD9);
    RETURN_CELL_TYPE_STR(t, HEX8);
    RETURN_CELL_TYPE_STR(t, HEX20);
    RETURN_CELL_TYPE_STR(t, HEX27);
    RETURN_CELL_TYPE_STR(t, TRI3);
    RETURN_CELL_TYPE_STR(t, TRI6);
    RETURN_CELL_TYPE_STR(t, TET4);
    RETURN_CELL_TYPE_STR(t, TET10);
    RETURN_CELL_TYPE_STR(t, PRISM6);
    RETURN_CELL_TYPE_STR(t, PRISM15);
    RETURN_CELL_TYPE_STR(t, PYRAMID5);
    RETURN_CELL_TYPE_STR(t, PYRAMID13);

    return "none";

#undef RETURN_CELL_TYPE_STR
}

const std::string MeshQualityType2String(const MeshQualityType t)
{
    if (t == MeshQualityType::ELEMENTSIZE)
        return "ElementSize";
    if (t == MeshQualityType::EDGERATIO)
        return "EdgeRatio";
    if (t == MeshQualityType::EQUIANGLESKEW)
        return "EquiAngleSkew";
    if (t == MeshQualityType::RADIUSEDGERATIO)
        return "RadiusEdgeRatio";
    if (t == MeshQualityType::SIZEDIFFERENCE)
        return "SizeDifference";
    return "none";
}

}
