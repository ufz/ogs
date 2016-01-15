/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_DIRICHLETBC_H
#define PROCESS_LIB_DIRICHLETBC_H

#include <vector>

#include "MeshLib/MeshSearch/NodeSearch.h"

namespace ProcessLib
{

/// A dirichlet boundary condition is represented by a list of global indices
/// with corresponding values.
template <typename IndexType>
struct DirichletBc final
{
	std::vector<IndexType> global_ids;
	std::vector<double> values;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_DIRICHLETBC_H
