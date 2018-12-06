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
            "The number of the given BHE properties (%d) are not consistent "
            "with the number of BHE groups in the mesh (%d).",
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
        auto const number_of_unknowns = apply_visitor(
            [](auto const& bhe) { return bhe.number_of_unknowns; },
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
    const double t, GlobalVector const& x, int const process_id)
{
    DBUG("Compute heat flux for HeatTransportBHE process.");
    GlobalExecutor::executeMemberOnDereferenced(
        &HeatTransportBHELocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, getDOFTable(process_id), t, x,
        _coupled_solutions);
}

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

        // Bottom and top nodes w.r.t. the z coordinate.
        auto const bottom_top_nodes = std::minmax_element(
            begin(bhe_nodes), end(bhe_nodes),
            [&](auto const& a, auto const& b) {
                return a->getCoords()[2] < b->getCoords()[2];
            });
        auto const bc_bottom_node_id = (*bottom_top_nodes.first)->getID();
        auto const bc_top_node_id = (*bottom_top_nodes.second)->getID();

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

        auto createBCs = [&](auto& bhe) {
            for (auto const& in_out_component_id :
                 bhe.inflow_outflow_bc_component_ids)
            {
                // Top, inflow.
                bcs.addBoundaryCondition(
                    ProcessLib::createBHEInflowDirichletBoundaryCondition(
                        get_global_bhe_bc_indices(bc_top_node_id,
                                                  in_out_component_id),
                        bhe));

                // Bottom, outflow.
                bcs.addBoundaryCondition(
                    ProcessLib::createBHEBottomDirichletBoundaryCondition(
                        get_global_bhe_bc_indices(bc_bottom_node_id,
                                                  in_out_component_id)));
            }
        };
        apply_visitor(createBCs, _process_data._vec_BHE_property[bhe_i]);
    }
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
