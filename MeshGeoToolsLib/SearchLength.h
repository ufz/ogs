/**
 * @file
 * @date 2014-09-19
 * @brief Base class for different search length strategies.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef SEARCHLENGTH_H_
#define SEARCHLENGTH_H_

#include "MeshLib/Mesh.h"

namespace MeshGeoToolsLib
{

/**
 * (Abstract) Base class for all SearchLength strategy implementations.
 */
class SearchLength
{
public:
	SearchLength(MeshLib::Mesh const& mesh)
		: _mesh(mesh), _search_length(0.0) {}
	virtual double getSearchLength() const = 0;
protected:
	MeshLib::Mesh const& _mesh;
	double _search_length;
};

} // end namespace MeshGeoToolsLib

#endif

