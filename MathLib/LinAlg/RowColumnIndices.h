/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

template <typename IDX_TYPE>
struct RowColumnIndices
{
	typedef typename std::vector<IDX_TYPE> LineIndex;
	RowColumnIndices(LineIndex const& rows, LineIndex const& columns)
		: rows(rows), columns(columns)
	{ }

	LineIndex const& rows;
	LineIndex const& columns;
};

} // MathLib

#endif  // ROWCOLUMNINDICES_H_
