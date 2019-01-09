/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEInflowPythonBoundaryCondition.h"

#include <pybind11/pybind11.h>
#include <iostream>

#include "MeshLib/MeshSearch/NodeSearch.h"
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
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
    int const variable_id,
    int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        pt_bhe,
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object)
    : _bc_mesh(bc_mesh), _pt_bhe(pt_bhe),_py_bc_object(py_bc_object)
{

    DBUG(
    "Found %d nodes for BHE Inflow Python BCs for the variable %d and "
    "component %d",
    vec_inflow_bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset{_bc_mesh, vec_inflow_bc_nodes};

    // create memory to store Tout values
    _T_out_values.clear();
    _T_out_indices.clear();

    _bc_values.ids.clear();
    _bc_values.values.clear();

        
    // convert mesh node ids to global index for the given component
    assert(bc_mesh_subset.getNumberOfNodes() == 1);
    _bc_values.ids.reserve(bc_mesh_subset.getNumberOfNodes());
    _bc_values.values.reserve(bc_mesh_subset.getNumberOfNodes());

    // that might be slow, but only done once
    const auto g_idx_T_in = in_out_global_indices.first;
    const auto g_idx_T_out = in_out_global_indices.second;

    if (g_idx_T_in >= 0 && g_idx_T_out >= 0)
    {
        _T_out_indices.emplace_back(g_idx_T_out);
        _T_out_values.emplace_back(320.0 /*using initial value*/);
        _bc_values.ids.emplace_back(g_idx_T_in);
        _bc_values.values.emplace_back(320.0 /*using initial value*/);
        //store the bc node ids to BHE network dataframe
        std::get<3>(_py_bc_object->dataframe_network).emplace_back(g_idx_T_out);
    }
}

void BHEInflowPythonBoundaryCondition::getEssentialBCValues(
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

void BHEInflowPythonBoundaryCondition::preTimestep(const double /*t*/,
                                                      const GlobalVector& x)
{
    // for each BHE, the inflow temperature is dependent on
    // the ouflow temperature of the BHE.
    // Here the task is to get the outflow temperature and
    // save it locally
    auto const n_nodes = _bc_values.ids.size();
    for (std::size_t i = 0; i < n_nodes; i++)
    {
        // read the T_out
        _T_out_values[i] = x[_T_out_indices[i]];
    }
}

void BHEInflowPythonBoundaryCondition::applyNaturalBC(const double t,
                                             const GlobalVector& x,
                                             GlobalMatrix& K, GlobalVector& b,
                                             GlobalMatrix* Jac)
{
}

std::unique_ptr<BHEInflowPythonBoundaryCondition>createBHEInflowPythonBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
    int const variable_id, int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const& pt_bhe,
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object)

{

    return std::make_unique<BHEInflowPythonBoundaryCondition>(
        std::move(in_out_global_indices), bc_mesh, vec_inflow_bc_nodes,
        variable_id, component_id, pt_bhe, py_bc_object);
}

}  // namespace ProcessLib
