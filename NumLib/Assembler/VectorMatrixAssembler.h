/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_VECTORMATRIXASSEMBLER_H_
#define NUMLIB_VECTORMATRIXASSEMBLER_H_

#include "NumLib/ODESolver/Types.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace NumLib
{

namespace detail
{
inline std::vector<GlobalIndexType>
getIndices(std::size_t const id,
           NumLib::LocalToGlobalIndexMap const& dof_table)
{
    assert(dof_table.size() > id);
    std::vector<GlobalIndexType> indices;

    // Local matrices and vectors will always be ordered by component
    // no matter what the order of the global matrix is.
    for (unsigned c = 0; c < dof_table.getNumberOfComponents(); ++c)
    {
        auto const& idcs = dof_table(id, c).rows;
        indices.reserve(indices.size() + idcs.size());
        indices.insert(indices.end(), idcs.begin(), idcs.end());
    }

    return indices;
}

inline std::vector<double>
getLocalNodalDOFs(GlobalVector const& x,
                  std::vector<GlobalIndexType> const& dof_indices)
{
    std::vector<double> local_x;
    local_x.reserve(dof_indices.size());

    for (auto i : dof_indices)
    {
        // TODO save some function calls to x[i]
        local_x.emplace_back(x[i]);
    }

    return local_x;
}
}  // namespace detail
}  // namespace NumLib

#endif  // NUMLIB_VECTORMATRIXASSEMBLER_H_
