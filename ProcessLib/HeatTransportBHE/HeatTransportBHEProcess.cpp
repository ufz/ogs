/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatTransportBHEProcess.h"

#include <cassert>

#include "BoundaryConditions/BHEBottomDirichletBoundaryCondition.h"
#include "BoundaryConditions/BHEInflowDirichletBoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/Python/BHEInflowPythonBoundaryCondition.h"
#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/CreateLocalAssemblers.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerBHE.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerSoil.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
HeatTransportBHEProcess::HeatTransportBHEProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HeatTransportBHEProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data)),
      _bheMeshData(getBHEDataInMesh(mesh))
{
    if (_bheMeshData.BHE_mat_IDs.size() !=
        _process_data._vec_BHE_property.size())
    {
        OGS_FATAL(
            "The number of the given BHE properties ({:d}) are not consistent "
            "with the number of BHE groups in the mesh ({:d}).",
            _process_data._vec_BHE_property.size(),
            _bheMeshData.BHE_mat_IDs.size());
    }

    auto material_ids = MeshLib::materialIDs(mesh);
    if (material_ids == nullptr)
    {
        OGS_FATAL("Not able to get material IDs! ");
    }

    _process_data._mesh_prop_materialIDs = material_ids;

    // create a map from a material ID to a BHE ID
    for (int i = 0; i < static_cast<int>(_bheMeshData.BHE_mat_IDs.size()); i++)
    {
        // fill in the map structure
        _process_data._map_materialID_to_BHE_ID[_bheMeshData.BHE_mat_IDs[i]] =
            i;
    }
}

void HeatTransportBHEProcess::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());

    //
    // Soil temperature variable defined on the whole mesh.
    //
    _mesh_subset_soil_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    std::vector<MeshLib::MeshSubset> all_mesh_subsets{*_mesh_subset_soil_nodes};

    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    vec_var_elements.push_back(&(_mesh.getElements()));

    std::vector<int> vec_n_components{
        1};  // one component for the soil temperature variable.

    //
    // BHE nodes with BHE type dependend number of variables.
    //
    int const n_BHEs = _process_data._vec_BHE_property.size();
    assert(n_BHEs == static_cast<int>(_bheMeshData.BHE_mat_IDs.size()));
    assert(n_BHEs == static_cast<int>(_bheMeshData.BHE_nodes.size()));
    assert(n_BHEs == static_cast<int>(_bheMeshData.BHE_elements.size()));

    // the BHE nodes need to be cherry-picked from the vector
    for (int i = 0; i < n_BHEs; i++)
    {
        auto const number_of_unknowns =
            visit([](auto const& bhe) { return bhe.number_of_unknowns; },
                  _process_data._vec_BHE_property[i]);
        auto const& bhe_nodes = _bheMeshData.BHE_nodes[i];
        auto const& bhe_elements = _bheMeshData.BHE_elements[i];

        // All the BHE nodes have additional variables.
        _mesh_subset_BHE_nodes.push_back(
            std::make_unique<MeshLib::MeshSubset const>(_mesh, bhe_nodes));

        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            // Here the number of components equals to the
            // number of unknowns on the BHE
            number_of_unknowns,
            [& ms = _mesh_subset_BHE_nodes.back()]() { return *ms; });

        vec_n_components.push_back(number_of_unknowns);
        vec_var_elements.push_back(&bhe_elements);
    }

    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    // in case of debugging the dof table, activate the following line
    // std::cout << *_local_to_global_index_map << "\n";
}

void HeatTransportBHEProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    // Quick access map to BHE's through element ids.
    std::unordered_map<std::size_t, BHE::BHETypes*> element_to_bhe_map;
    int const n_BHEs = _process_data._vec_BHE_property.size();
    for (int i = 0; i < n_BHEs; i++)
    {
        auto const& bhe_elements = _bheMeshData.BHE_elements[i];
        for (auto const& e : bhe_elements)
        {
            element_to_bhe_map[e->getID()] =
                &_process_data._vec_BHE_property[i];
        }
    }

    assert(mesh.getDimension() == 3);
    ProcessLib::HeatTransportBHE::createLocalAssemblers<
        HeatTransportBHELocalAssemblerSoil, HeatTransportBHELocalAssemblerBHE>(
        mesh.getElements(), dof_table, _local_assemblers, element_to_bhe_map,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    // Create BHE boundary conditions for each of the BHEs
    createBHEBoundaryConditionTopBottom(_bheMeshData.BHE_nodes);
}

void HeatTransportBHEProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble HeatTransportBHE process.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b);
}

void HeatTransportBHEProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
    const double dx_dx, int const process_id, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HeatTransportBHE process.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac);
}

void HeatTransportBHEProcess::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_dot, int const process_id)
{
    DBUG("Compute heat flux for HeatTransportBHE process.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HeatTransportBHELocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_dot, process_id);
}

#ifdef OGS_USE_PYTHON
NumLib::IterationResult HeatTransportBHEProcess::postIterationConcreteProcess(
    GlobalVector const& x)
{
    // if the process use python boundary conditon
    if (_process_data.py_bc_object == nullptr)
        return NumLib::IterationResult::SUCCESS;

    // Here the task is to get current time flowrate and flow temperature from
    // TESPy and determine whether it converges.
    auto const Tout_nodes_id =
        std::get<3>(_process_data.py_bc_object->dataframe_network);
    const std::size_t n_bc_nodes = Tout_nodes_id.size();

    for (std::size_t i = 0; i < n_bc_nodes; i++)
    {
        // read the T_out and store them in dataframe
        std::get<2>(_process_data.py_bc_object->dataframe_network)[i] =
            x[Tout_nodes_id[i]];
    }
    // Transfer Tin and Tout to TESPy and return the results
    auto const tespy_result = _process_data.py_bc_object->tespySolver(
        std::get<0>(_process_data.py_bc_object->dataframe_network),
        std::get<1>(_process_data.py_bc_object->dataframe_network),
        std::get<2>(_process_data.py_bc_object->dataframe_network));
    if (!_process_data.py_bc_object->isOverriddenTespy())
    {
        DBUG(
            "Method `tespySolver' not overridden in Python "
            "script.");
    }

    // update the Tin and flow rate
    for (std::size_t i = 0; i < n_bc_nodes; i++)
    {
        std::get<1>(_process_data.py_bc_object->dataframe_network)[i] =
            std::get<2>(tespy_result)[i];
        std::get<4>(_process_data.py_bc_object->dataframe_network)[i] =
            std::get<3>(tespy_result)[i];
    }
    auto const tespy_has_converged = std::get<1>(tespy_result);
    if (tespy_has_converged == true)
        return NumLib::IterationResult::SUCCESS;

    return NumLib::IterationResult::REPEAT_ITERATION;
}
#endif  // OGS_USE_PYTHON

void HeatTransportBHEProcess::createBHEBoundaryConditionTopBottom(
    std::vector<std::vector<MeshLib::Node*>> const& all_bhe_nodes)
{
    const int process_id = 0;
    auto& bcs = _boundary_conditions[process_id];

    int const n_BHEs = static_cast<int>(_process_data._vec_BHE_property.size());

    // for each BHE
    for (int bhe_i = 0; bhe_i < n_BHEs; bhe_i++)
    {
        auto const& bhe_nodes = all_bhe_nodes[bhe_i];
        // find the variable ID
        // the soil temperature is 0-th variable
        // the BHE temperature is therefore bhe_i + 1
        const int variable_id = bhe_i + 1;

        std::vector<MeshLib::Node*> bhe_boundary_nodes;

        // cherry-pick the boundary nodes according to
        // the number of connected line elements.
        for (auto const& bhe_node : bhe_nodes)
        {
            // Count number of 1d elements connected with every BHE node.
            const std::size_t n_line_elements = std::count_if(
                bhe_node->getElements().begin(), bhe_node->getElements().end(),
                [](MeshLib::Element const* elem) {
                    return (elem->getDimension() == 1);
                });

            if (n_line_elements == 1)
            {
                bhe_boundary_nodes.push_back(bhe_node);
            }
        }

        if (bhe_boundary_nodes.size() != 2)
        {
            OGS_FATAL(
                "Error!!! The BHE boundary nodes are not correctly found, "
                "for every single BHE, there should be 2 boundary nodes.");
        }

        // For 1U, 2U, CXC, CXA type BHE, the node order in the boundary nodes
        // vector should be rearranged according to its z coordinate in
        // descending order. In these BHE types, the z coordinate on the top and
        // bottom node is different. The BHE top node with a higher z coordinate
        // should be placed at the first, while the BHE bottom node with a lower
        // z coordinate should be placed at the second. For other horizontal BHE
        // types e.g. 1P-type BHE, the z coordinate on the top and bottom node
        // is identical. Thus the node order in the boundary nodes vector can
        // not be rearranged according to its z coordinate. For these BHE types,
        // the boundary node order is according to the default node id order in
        // the model mesh.
        // for 1P-type BHE
        if (bhe_boundary_nodes[0]->getCoords()[2] ==
            bhe_boundary_nodes[1]->getCoords()[2])
        {
            INFO(
                "For 1P-type BHE, the BHE inflow and outflow "
                "nodes are identified according to their mesh node "
                "id in ascending order");
        }
        // for 1U, 2U, CXC, CXA type BHE
        else
        {
            // swap the boundary nodes if the z coordinate of the
            // first node is lower than it on the second node
            if (bhe_boundary_nodes[0]->getCoords()[2] <
                bhe_boundary_nodes[1]->getCoords()[2])
            {
                std::swap(bhe_boundary_nodes[0], bhe_boundary_nodes[1]);
            }
        }

        auto get_global_index =
            [&](std::size_t const node_id, int const component) {
                return _local_to_global_index_map->getGlobalIndex(
                    {_mesh.getID(), MeshLib::MeshItemType::Node, node_id},
                    variable_id, component);
            };

        auto get_global_bhe_bc_indices =
            [&](std::array<
                std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
                    nodes_and_components) {
                return std::make_pair(
                    get_global_index(nodes_and_components[0].first,
                                     nodes_and_components[0].second),
                    get_global_index(nodes_and_components[1].first,
                                     nodes_and_components[1].second));
            };

        auto createBCs = [&, bc_top_node_id = bhe_boundary_nodes[0]->getID(),
                          bc_bottom_node_id =
                              bhe_boundary_nodes[1]->getID()](auto& bhe) {
            for (auto const& in_out_component_id :
                 bhe.inflow_outflow_bc_component_ids)
            {
                if (bhe.use_python_bcs)
                // call BHEPythonBoundarycondition
                {
#ifdef OGS_USE_PYTHON
                    if (_process_data.py_bc_object)  // the bc object exist
                    {
                        // apply the customized top, inflow BC.
                        bcs.addBoundaryCondition(
                            ProcessLib::createBHEInflowPythonBoundaryCondition(
                                get_global_bhe_bc_indices(
                                    bhe.getBHEInflowDirichletBCNodesAndComponents(
                                        bc_top_node_id, bc_bottom_node_id,
                                        in_out_component_id.first)),
                                bhe,
                                *(_process_data.py_bc_object)));
                    }
                    else
                    {
                        OGS_FATAL(
                            "The Python Boundary Condition was switched on, "
                            "but the data object does not exist! ");
                    }
#else
                    OGS_FATAL(
                        "The Python Boundary Condition was switched off! "
                        "Not able to create Boundary Condtion for BHE! ");
#endif  // OGS_USE_PYTHON
                }
                else
                {
                    // Top, inflow, normal case
                    bcs.addBoundaryCondition(
                        createBHEInflowDirichletBoundaryCondition(
                            get_global_bhe_bc_indices(
                                bhe.getBHEInflowDirichletBCNodesAndComponents(
                                    bc_top_node_id, bc_bottom_node_id,
                                    in_out_component_id.first)),
                            [&bhe](double const T, double const t) {
                                return bhe.updateFlowRateAndTemperature(T, t);
                            }));
                }

                auto const bottom_nodes_and_components =
                    bhe.getBHEBottomDirichletBCNodesAndComponents(
                        bc_bottom_node_id,
                        in_out_component_id.first,
                        in_out_component_id.second);

                if (bottom_nodes_and_components)
                {
                    // Bottom, outflow, all cases
                    bcs.addBoundaryCondition(
                        createBHEBottomDirichletBoundaryCondition(
                            get_global_bhe_bc_indices(
                                {{{bc_bottom_node_id, in_out_component_id.first},
                                  {bc_bottom_node_id,
                                   in_out_component_id.second}}})));
                }

            }
        };
        visit(createBCs, _process_data._vec_BHE_property[bhe_i]);
    }
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
