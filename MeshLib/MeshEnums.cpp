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
	if (t == MeshElemType::EDGE)
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
	if (s.compare("line") == 0)
		return MeshElemType::EDGE;
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

const std::string MshQualityType2String(const MshQualityType::type t)
{
	if (t == MshQualityType::AREA)
		return "Area";
	if (t == MshQualityType::EDGERATIO)
		return "EdgeRatio";
	if (t == MshQualityType::EQUIANGLESKEW)
		return "EquiAngleSkew";
	if (t == MshQualityType::VOLUME)
		return "Volume";
	return "none";
}
