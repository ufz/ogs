/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of mesh-related enumerations.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshEnums.h"

const std::string MeshElemType2String(const MeshElemType t)
{
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

MeshElemType CompleteString2MeshElemType(const std::string &s)
{
	if (s.compare("Line") == 0)
		return MeshElemType::LINE;
	if (s.compare("Quad") == 0)
		return MeshElemType::QUAD;
	if (s.compare("Hexahedron")  == 0)
		return MeshElemType::HEXAHEDRON;
	if (s.compare("Triangle")  == 0)
		return MeshElemType::TRIANGLE;
	if (s.compare("Tetrahedron")  == 0)
		return MeshElemType::TETRAHEDRON;
	if (s.compare("Prism") == 0)
		return MeshElemType::PRISM;
	if (s.compare("Pyramid") == 0)
		return MeshElemType::PYRAMID;
	return MeshElemType::INVALID;
}

MeshElemType String2MeshElemType(const std::string &s)
{
	if (s.compare("line") == 0)
		return MeshElemType::LINE;
	if (s.compare("quad") == 0)
		return MeshElemType::QUAD;
	if (s.compare("hex")  == 0)
		return MeshElemType::HEXAHEDRON;
	if (s.compare("tri")  == 0)
		return MeshElemType::TRIANGLE;
	if (s.compare("tet")  == 0)
		return MeshElemType::TETRAHEDRON;
	if (s.compare("pris") == 0)
		return MeshElemType::PRISM;
	if (s.compare("pyra") == 0)
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
