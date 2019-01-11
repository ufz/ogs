/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEInflowPythonBoundaryCondition.h"
#include "BaseLib/Error.h"
#include <pybind11/pybind11.h>
#include <iostream>

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "PythonBoundaryConditionLocalAssembler.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>



namespace ProcessLib
{
BHEInflowPythonBoundaryCondition::BHEInflowPythonBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object)
    : _in_out_global_indices(std::move(in_out_global_indices)),
      _py_bc_object(py_bc_object)
{

    const auto g_idx_T_out = in_out_global_indices.second;

    //store the bc node ids to BHE network dataframe
    std::get<3>(_py_bc_object->dataframe_network).emplace_back(g_idx_T_out);
    
}

void BHEInflowPythonBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    bc_values.ids.resize(1);
    bc_values.values.resize(1);

    bc_values.ids[0] = _in_out_global_indices.first;
    // here call the corresponding BHE functions
    auto const T_out = x[_in_out_global_indices.second];

}

std::unique_ptr<BHEInflowPythonBoundaryCondition>
createBHEInflowPythonBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object)

{
    DBUG("Constructing BHEInflowPythonBoundaryCondition.");

    // In case of partitioned mesh the boundary could be empty, i.e. there is no
    // boundary condition.
#ifdef USE_PETSC
    // For this special boundary condition the boundary condition is not empty
    // if the global indices are non-negative.
    if (in_out_global_indices.first < 0 && in_out_global_indices.second < 0)
    {
        return nullptr;
    }
    // If only one of the global indices (in or out) is negative the
    // implementation is not valid.
    if (in_out_global_indices.first < 0 || in_out_global_indices.second < 0)
    {
        OGS_FATAL(
            "The partition cuts the BHE into two independent parts. This "
            "behaviour is not implemented.");
    }
#endif  // USE_PETSC
    return std::make_unique<BHEInflowPythonBoundaryCondition>(
        std::move(in_out_global_indices),  py_bc_object);
}

}  // namespace ProcessLib
