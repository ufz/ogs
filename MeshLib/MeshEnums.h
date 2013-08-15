/**
 * \file
 * \author  Karsten Rink
 * \date    2012-05-02
 * \brief   Definition of mesh-related Enumerations.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHENUMS_H
#define MESHENUMS_H

#include <string>

/**
 * \brief Types of mesh elements supported by OpenGeoSys.
 */
enum class MeshElemType
{
	INVALID,
	EDGE,
	QUAD,
	HEXAHEDRON,
	TRIANGLE,
	TETRAHEDRON,
	PRISM,
	PYRAMID
};

/**
 * \brief Types of mesh elements supported by OpenGeoSys.
 */
struct CellType
{
	enum type {
		INVALID,
		EDGE2,
		EDGE3,
		TRI3,
		TRI6,
		QUAD4,
		QUAD8,
		QUAD9,
		TET4,
		TET10,
		HEX8,
		HEX20,
		HEX27,
		PRISM6,
		PRISM15,
		PRISM18,
		PYRAMID5
	};
};

/**
 * \brief Describes a mesh quality criteria.
 */
struct MshQualityType
{
	enum type {
		INVALID = 0,
		AREA,
		VOLUME,
		EDGERATIO,
		EQUIANGLESKEW
	};
};

/// Given a MeshElemType this returns the appropriate string.
const std::string MeshElemType2String(const MeshElemType t);

/// Given a string describing an element type this returns the corresponding MeshElemType.
MeshElemType String2MeshElemType(const std::string &s);

const std::string MshQualityType2String(const MshQualityType::type t);

#endif //MESHENUMS_H
