/**
 * @file
 * @date 2014-09-22
 * @brief Interface / impl. for fixed search length strategy.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef FIXEDSEARCHLENGTH_H_
#define FIXEDSEARCHLENGTH_H_

#include <limits>

#include "MeshGeoToolsLib/SearchLength.h"

namespace MeshGeoToolsLib
{

/// FixedSearchLength implements a mesh independent, strong criterion for
/// searching mesh nodes near a geometry. The algorithm can be used for meshes
/// that have nearly equi-sized elements.
class FixedSearchLength : public SearchLength
{
public:
	/// Constructor for FixedSearchLength object with a default search length
	/// of 10 angstrom (\f$10^{-9}\f$ m)
	FixedSearchLength(MeshLib::Mesh const& mesh,
		double search_length = 1e-9)
	: SearchLength(mesh, search_length)
	{}

	virtual double getSearchLength() const { return _search_length; }
};

} // end namespace MeshGeoToolsLib

#endif

