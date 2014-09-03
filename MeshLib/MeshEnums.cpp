/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of mesh-related enumerations.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshEnums.h"

const std::string MeshElemType2String(const MeshElemType t)
{
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

MeshElemType String2MeshElemType(const std::string &s)
{
	if ((s.compare("line") == 0) || (s.compare("Line") == 0))
		return MeshElemType::LINE;
	if ((s.compare("quad") == 0) || (s.compare("Quadrilateral") == 0))
		return MeshElemType::QUAD;
	if ((s.compare("hex")  == 0) || (s.compare("Hexahedron") == 0))
		return MeshElemType::HEXAHEDRON;
	if ((s.compare("tri")  == 0) || (s.compare("Triangle") == 0))
		return MeshElemType::TRIANGLE;
	if ((s.compare("tet")  == 0) || (s.compare("Tetrahedron") == 0))
		return MeshElemType::TETRAHEDRON;
	if ((s.compare("pris") == 0) || (s.compare("Prism") == 0))
		return MeshElemType::PRISM;
	if ((s.compare("pyra") == 0) || (s.compare("Pyramid") == 0))
		return MeshElemType::PYRAMID;
	return MeshElemType::INVALID;
}

const std::string MeshQualityType2String(const MeshQualityType t)
{
	if (t == MeshQualityType::AREA)
		return "Area";
	if (t == MeshQualityType::EDGERATIO)
		return "EdgeRatio";
	if (t == MeshQualityType::EQUIANGLESKEW)
		return "EquiAngleSkew";
	if (t == MeshQualityType::VOLUME)
		return "Volume";
	return "none";
}

