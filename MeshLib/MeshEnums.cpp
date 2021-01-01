/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of mesh-related enumerations.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshEnums.h"

#include <boost/algorithm/string/predicate.hpp>

namespace MeshLib {
std::string MeshElemType2String(const MeshElemType t)
{
    if (t == MeshElemType::POINT)
    {
        return "Point";
    }
    if (t == MeshElemType::LINE)
    {
        return "Line";
    }
    if (t == MeshElemType::QUAD)
    {
        return "Quad";
    }
    if (t == MeshElemType::HEXAHEDRON)
    {
        return "Hexahedron";
    }
    if (t == MeshElemType::TRIANGLE)
    {
        return "Triangle";
    }
    if (t == MeshElemType::TETRAHEDRON)
    {
        return "Tetrahedron";
    }
    if (t == MeshElemType::PRISM)
    {
        return "Prism";
    }
    if (t == MeshElemType::PYRAMID)
    {
        return "Pyramid";
    }
    return "none";
}

std::string MeshElemType2StringShort(const MeshElemType t)
{
    if (t == MeshElemType::POINT)
    {
        return "point";
    }
    if (t == MeshElemType::LINE)
    {
        return "line";
    }
    if (t == MeshElemType::QUAD)
    {
        return "quad";
    }
    if (t == MeshElemType::HEXAHEDRON)
    {
        return "hex";
    }
    if (t == MeshElemType::TRIANGLE)
    {
        return "tri";
    }
    if (t == MeshElemType::TETRAHEDRON)
    {
        return "tet";
    }
    if (t == MeshElemType::PRISM)
    {
        return "pris";
    }
    if (t == MeshElemType::PYRAMID)
    {
        return "pyra";
    }
    return "none";
}

MeshElemType String2MeshElemType(const std::string &s)
{
    if (boost::iequals(s, "point"))
    {
        return MeshElemType::POINT;
    }
    if (boost::iequals(s, "line"))
    {
        return MeshElemType::LINE;
    }
    if (boost::iequals(s, "quad") || boost::iequals(s, "Quadrilateral"))
    {
        return MeshElemType::QUAD;
    }
    if (boost::iequals(s, "hex") || boost::iequals(s, "Hexahedron"))
    {
        return MeshElemType::HEXAHEDRON;
    }
    if (boost::iequals(s, "tri") || boost::iequals(s, "Triangle"))
    {
        return MeshElemType::TRIANGLE;
    }
    if (boost::iequals(s, "tet") || boost::iequals(s, "Tetrahedron"))
    {
        return MeshElemType::TETRAHEDRON;
    }
    if (boost::iequals(s, "pris") || boost::iequals(s, "Prism"))
    {
        return MeshElemType::PRISM;
    }
    if (boost::iequals(s, "pyra") || boost::iequals(s, "Pyramid"))
    {
        return MeshElemType::PYRAMID;
    }
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
    auto const& mesh_elem_types = getMeshElemTypes();
    std::transform(mesh_elem_types.begin(), mesh_elem_types.end(),
                   std::back_inserter(vec), [](auto const& element_type) {
                       return MeshElemType2StringShort(element_type);
                   });
    return vec;
}

std::string CellType2String(const CellType t)
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

std::string MeshQualityType2String(const MeshQualityType t)
{
    if (t == MeshQualityType::ELEMENTSIZE)
    {
        return "ElementSize";
    }
    if (t == MeshQualityType::EDGERATIO)
    {
        return "EdgeRatio";
    }
    if (t == MeshQualityType::EQUIANGLESKEW)
    {
        return "EquiAngleSkew";
    }
    if (t == MeshQualityType::RADIUSEDGERATIO)
    {
        return "RadiusEdgeRatio";
    }
    if (t == MeshQualityType::SIZEDIFFERENCE)
    {
        return "SizeDifference";
    }
    return "none";
}

MeshLib::MeshQualityType String2MeshQualityType(std::string const& s)
{
    if (boost::iequals(s, "ElementSize"))
    {
        return MeshQualityType::ELEMENTSIZE;
    }
    if (boost::iequals(s, "EdgeRatio"))
    {
        return MeshQualityType::EDGERATIO;
    }
    if (boost::iequals(s, "EquiAngleSkew"))
    {
        return MeshQualityType::EQUIANGLESKEW;
    }
    if (boost::iequals(s, "RadiusEdgeRatio"))
    {
        return MeshQualityType::RADIUSEDGERATIO;
    }
    if (boost::iequals(s, "SizeDifference"))
    {
        return MeshQualityType::SIZEDIFFERENCE;
    }
    return MeshQualityType::INVALID;
}

}  // namespace MeshLib
