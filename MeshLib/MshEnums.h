/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
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
