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

const std::string MshElemType2String(const MshElemType::type t)
{
	if (t == MshElemType::EDGE)
		return "line";
	if (t == MshElemType::QUAD)
		return "quad";
	if (t == MshElemType::HEXAHEDRON)
		return "hex";
	if (t == MshElemType::TRIANGLE)
		return "tri";
	if (t == MshElemType::TETRAHEDRON)
		return "tet";
	if (t == MshElemType::PRISM)
		return "pris";
	if (t == MshElemType::PYRAMID)
		return "pyra";
	return "none";
}

MshElemType::type String2MshElemType(const std::string &s)
{
	if (s.compare("line") == 0)
		return MshElemType::EDGE;
	if (s.compare("quad") == 0)
		return MshElemType::QUAD;
	if (s.compare("hex")  == 0)
		return MshElemType::HEXAHEDRON;
	if (s.compare("tri")  == 0)
		return MshElemType::TRIANGLE;
	if (s.compare("tet")  == 0)
		return MshElemType::TETRAHEDRON;
	if (s.compare("pris") == 0)
		return MshElemType::PRISM;
	if (s.compare("pyra") == 0)
		return MshElemType::PYRAMID;
	return MshElemType::INVALID;
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
