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

#include "PythonBoundaryConditionLocalAssembler.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>



namespace ProcessLib
{
template <typename BHEType>
BHEInflowPythonBoundaryCondition<BHEType>::BHEInflowPythonBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    BHEType& bhe,
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object)
    : _in_out_global_indices(std::move(in_out_global_indices)),
      _bhe(bhe),
      _py_bc_object(py_bc_object)
{

    const auto g_idx_T_out = in_out_global_indices.second;

    //store the bc node ids to BHE network dataframe
    std::get<3>(_py_bc_object->dataframe_network).emplace_back(g_idx_T_out);
    
}

template <typename BHEType>
void BHEInflowPythonBoundaryCondition<BHEType>::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.resize(_bc_values.ids.size());
    bc_values.values.resize(_bc_values.values.size());

    const std::size_t n_nodes = _T_out_values.size();

    std::vector<double> primary_variables;

    // get the number of all boundary nodes
    const std::size_t n_bc_nodes = std::get<3>(_py_bc_object->dataframe_network).size();

    for (std::size_t i = 0; i < n_nodes; i++)
    {
        bc_values.ids[i] = _bc_values.ids[i];

        // here call the corresponding BHE functions
        auto const tmp_T_out = x[_T_out_indices[i]];

        primary_variables.push_back(tmp_T_out);

        auto const boundary_node_id = _T_out_indices[i];

        //return T_in from currently BHE dataframe column 2
        for (std::size_t j = 0; j < n_bc_nodes; j++)
        {
            auto const dataframe_node_id = std::get<3>(_py_bc_object->dataframe_network);
            auto const dataframe_Tin_val = std::get<1>(_py_bc_object->dataframe_network);
            if (dataframe_node_id[j] == boundary_node_id)
            {
                bc_values.values[i] = dataframe_Tin_val[j];
                break;
            }
        }
    }

}

template <typename BHEType>
std::unique_ptr<BHEInflowPythonBoundaryCondition<BHEType>>
createBHEInflowPythonBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    BHEType& bhe,
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
    return std::make_unique<BHEInflowPythonBoundaryCondition<BHEType>>(
        std::move(in_out_global_indices), bhe, py_bc_object);
}

}  // namespace ProcessLib
