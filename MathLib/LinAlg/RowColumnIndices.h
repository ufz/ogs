/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

#ifdef USE_PETSC
template <typename IDX_TYPE>
struct RowColumnIndices
{
	typedef typename std::vector<IDX_TYPE> LineIndex;
	RowColumnIndices(LineIndex const& rows_, LineIndex const& columns_, const std::vector<bool>& ghost_flags_)
		: rows(rows_), columns(columns_), ghost_flags(ghost_flags_)
	{ }

	LineIndex const& rows;
	LineIndex const& columns;
   
    std::vector<bool> const& ghost_flags;
};
#else
template <typename IDX_TYPE>
struct RowColumnIndices
{
	typedef typename std::vector<IDX_TYPE> LineIndex;
	RowColumnIndices(LineIndex const& rows_, LineIndex const& columns_)
		: rows(rows_), columns(columns_)
	{ }

	LineIndex const& rows;
	LineIndex const& columns;
};
#endif
} // MathLib

#endif  // ROWCOLUMNINDICES_H_
