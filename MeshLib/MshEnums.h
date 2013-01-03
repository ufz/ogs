/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file MshEnums.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef MSHENUMS_H
#define MSHENUMS_H

#include <string>

/**
 * \brief Types of mesh elements supported by OpenGeoSys.
 */
struct MshElemType
{
	enum type {
		INVALID,
		EDGE,
		QUAD,
		HEXAHEDRON,
		TRIANGLE,
		TETRAHEDRON,
		PRISM,
		PYRAMID
	};
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

/// Given a MshElemType this returns the appropriate string.
const std::string MshElemType2String(const MshElemType::type t);

/// Given a string describing an element type this returns the corresponding MshElemType.
MshElemType::type String2MshElemType(const std::string &s);

const std::string MshQualityType2String(const MshQualityType::type t);

#endif //MSHENUMS_H
