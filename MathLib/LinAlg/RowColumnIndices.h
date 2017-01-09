/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    RowColumnIndices(LineIndex const& rows_, LineIndex const& columns_)
        : rows(rows_), columns(columns_)
    { }

    LineIndex const& rows;
    LineIndex const& columns;
};

} // MathLib

#endif  // ROWCOLUMNINDICES_H_
