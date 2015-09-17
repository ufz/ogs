/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ROWCOLUMNINDICES_H_
#define ROWCOLUMNINDICES_H_

#include <vector>

namespace MathLib
{
// A dummy vector as the default value of the third argument of
// the constructor.
static std::vector<std::size_t> dummy_ghost_local_ids;

template <typename IDX_TYPE>
struct RowColumnIndices
{
	typedef typename std::vector<IDX_TYPE> LineIndex;
	RowColumnIndices(LineIndex const& rows_, LineIndex const& columns_,
	                 const std::vector<std::size_t>& ghost_local_ids_
	                 = dummy_ghost_local_ids)
		: rows(rows_), columns(columns_), ghost_local_ids(ghost_local_ids_)
	{ }

	LineIndex const& rows;
	LineIndex const& columns;

	std::vector<std::size_t> const& ghost_local_ids;
};
} // MathLib

#endif  // ROWCOLUMNINDICES_H_
