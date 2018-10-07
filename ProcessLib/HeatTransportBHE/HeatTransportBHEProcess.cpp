/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatTransportBHEProcess.h"

#include <cassert>
#include <iostream>
#include "ProcessLib/HeatTransportBHE/BHE/MeshUtils.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/CreateLocalAssemblers.h"

#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerBHE.h"
#include "ProcessLib/HeatTransportBHE/LocalAssemblers/HeatTransportBHELocalAssemblerSoil.h"

#include "ProcessLib/BoundaryCondition/BHEBottomDirichletBoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/BHEInflowDirichletBoundaryCondition.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
HeatTransportBHEProcess::HeatTransportBHEProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HeatTransportBHEProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data)),
      _bheMeshData(getBHEDataInMesh(mesh))
{
    if (_bheMeshData.BHE_mat_IDs.size() !=
        _process_data._vec_BHE_property.size())
    {
        OGS_FATAL(
            "The number of the given BHE properties (%d) are not "
            "consistent"
            " with the number of BHE groups in a mesh (%d).",
            _process_data._vec_BHE_property.size(),
            _bheMeshData.BHE_mat_IDs.size());
    }

    auto material_ids = MeshLib::materialIDs(mesh);
    if (material_ids == nullptr)
    {
        OGS_FATAL("Not able to get material IDs! ");
    }
    // create a map from a material ID to a BHE ID
    auto max_BHE_mat_id =
        std::max_element(material_ids->begin(), material_ids->end());
    _process_data._map_materialID_to_BHE_ID.resize(*max_BHE_mat_id + 1);
    for (std::size_t i = 0; i < _bheMeshData.BHE_mat_IDs.size(); i++)
    {
        // fill in the map structure
        _process_data._map_materialID_to_BHE_ID[_bheMeshData.BHE_mat_IDs[i]] =
            i;
    }

    _process_data._mesh_prop_materialIDs = material_ids;
}

void HeatTransportBHEProcess::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());

    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // first all the soil nodes
    _mesh_subset_pure_soil_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _bheMeshData.soil_nodes);

    std::vector<int> vec_n_BHE_unknowns;
    vec_n_BHE_unknowns.reserve(_bheMeshData.BHE_nodes.size());
    // the BHE nodes need to be cherry-picked from the vector
    for (unsigned i = 0; i < _bheMeshData.BHE_nodes.size(); i++)
    {
        _mesh_subset_BHE_nodes.push_back(
            std::make_unique<MeshLib::MeshSubset const>(
                _mesh, _bheMeshData.BHE_nodes[i]));
        int n_BHE_unknowns =
            _process_data._vec_BHE_property[i]->getNumUnknowns();
        vec_n_BHE_unknowns.emplace_back(n_BHE_unknowns);
    }

    // Collect the mesh subsets in a vector.
    // All the soil nodes has 1 temperature variable
    std::vector<MeshLib::MeshSubset> all_mesh_subsets{
        *_mesh_subset_pure_soil_nodes};

    // All the BHE nodes have additinal variables
    std::size_t count = 0;
    for (auto& ms : _mesh_subset_BHE_nodes)
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        // Here the number of components equals to
                        // the number of unknowns on the BHE
                        vec_n_BHE_unknowns[count],
                        [&]() { return *ms; });
        count++;
    }

    std::vector<int> vec_n_components;
    // This is the soil temperature for first mesh subset.
    // 1, because for the soil part there is just one variable which is the soil
    // temperature.
    vec_n_components.push_back(1);
    // now the BHE subsets
    for (std::size_t i = 0; i < _bheMeshData.BHE_mat_IDs.size(); i++)
    {
        // Here the number of components equals to
        // the number of unknowns on the BHE
        vec_n_components.push_back(vec_n_BHE_unknowns[i]);
    }

    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    // vec_var_elements.push_back(&_vec_pure_soil_elements);
    vec_var_elements.push_back(&(_mesh.getElements()));
    for (std::size_t i = 0; i < _bheMeshData.BHE_elements.size(); i++)
    {
        vec_var_elements.push_back(&_bheMeshData.BHE_elements[i]);
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
    assert(mesh.getDimension() == 3);
    ProcessLib::HeatTransportBHE::createLocalAssemblers<
        3 /* mesh dimension */, HeatTransportBHELocalAssemblerSoil,
        HeatTransportBHELocalAssemblerBHE>(
        mesh.getElements(), dof_table, _local_assemblers,
        _process_data._vec_ele_connected_BHE_IDs,
        _process_data._vec_BHE_property, mesh.isAxiallySymmetric(),
        integration_order, _process_data);

    // create BHE boundary conditions
    // for each BHE, one BC on the top
    // and one BC at the bottom.
    std::vector<std::unique_ptr<BoundaryCondition>> bc_collections =
        createBHEBoundaryConditionTopBottom();

    const int process_id = 0;
    auto& current_process_BCs = _boundary_conditions[process_id];
    for (auto& bc_collection : bc_collections)
    {
        current_process_BCs.addCreatedBC(std::move(bc_collection));
    }
}

void HeatTransportBHEProcess::assembleConcreteProcess(const double t,
                                                      GlobalVector const& x,
                                                      GlobalMatrix& M,
                                                      GlobalMatrix& K,
                                                      GlobalVector& b)
{
    DBUG("Assemble HeatTransportBHE process.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

void HeatTransportBHEProcess::assembleWithJacobianConcreteProcess(
    const double /*t*/, GlobalVector const& /*x*/, GlobalVector const& /*xdot*/,
    const double /*dxdot_dx*/, const double /*dx_dx*/, GlobalMatrix& /*M*/,
    GlobalMatrix& /*K*/, GlobalVector& /*b*/, GlobalMatrix& /*Jac*/)
{
    OGS_FATAL(
        "HeatTransportBHE: analytical Jacobian assembly is not implemented");
}

void HeatTransportBHEProcess::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x)
{
    DBUG("Compute heat flux for HeatTransportBHE process.");
    GlobalExecutor::executeMemberOnDereferenced(
        &HeatTransportBHELocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, *_local_to_global_index_map, t, x,
        _coupled_solutions);
}

std::vector<std::unique_ptr<BoundaryCondition>>
HeatTransportBHEProcess::createBHEBoundaryConditionTopBottom()
{
    /**
     * boundary conditions
     */
    std::vector<std::unique_ptr<BoundaryCondition>> bcs;

    int const n_BHEs = static_cast<int>(_process_data._vec_BHE_property.size());

    // for each BHE
    for (int bhe_i = 0; bhe_i < n_BHEs; bhe_i++)
    {
        auto const& bhe_nodes = _bheMeshData.BHE_nodes[bhe_i];
        // find the variable ID
        // the soil temperature is 0-th variable
        // the BHE temperature is therefore bhe_i + 1
        const int variable_id = bhe_i + 1;

        // Bottom and top nodes w.r.t. the z coordinate.
        auto const bottom_top_nodes = std::minmax_element(
            begin(bhe_nodes), end(bhe_nodes),
            [&](auto const& a, auto const& b) {
                return a->getCoords()[2] < b->getCoords()[2];
            });
        auto const bc_bottom_node_id = (*bottom_top_nodes.first)->getID();
        MeshLib::Node* const bc_top_node = *bottom_top_nodes.second;

        auto get_global_bhe_bc_indices =
            [&](std::size_t const node_id,
                std::pair<int, int> const& in_out_component_id) {
                return std::make_pair(
                    _local_to_global_index_map->getGlobalIndex(
                        {_mesh.getID(), MeshLib::MeshItemType::Node, node_id},
                        variable_id, in_out_component_id.first),
                    _local_to_global_index_map->getGlobalIndex(
                        {_mesh.getID(), MeshLib::MeshItemType::Node, node_id},
                        variable_id, in_out_component_id.second));
            };

        for (auto const& in_out_component_id :
             _process_data._vec_BHE_property[bhe_i]
                 ->inflowOutflowBcComponentIds())
        {
            // Top, inflow.
            bcs.push_back(ProcessLib::createBHEInflowDirichletBoundaryCondition(
                get_global_bhe_bc_indices(bc_top_node->getID(),
                                          in_out_component_id),
                _mesh, {bc_top_node}, variable_id, in_out_component_id.first,
                _process_data._vec_BHE_property[bhe_i]));

            // Bottom, outflow.
            bcs.push_back(ProcessLib::createBHEBottomDirichletBoundaryCondition(
                get_global_bhe_bc_indices(bc_bottom_node_id,
                                          in_out_component_id)));
        }
    }

    return bcs;
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
