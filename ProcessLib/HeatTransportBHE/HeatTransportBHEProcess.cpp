/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HeatTransportBHEProcess.h"

#include <cassert>

#include "BoundaryConditions/BHEBottomDirichletBoundaryCondition.h"
#include "BoundaryConditions/BHEInflowDirichletBoundaryCondition.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/BHEInflowPythonBoundaryCondition.h"
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
    SecondaryVariableCollection&& secondary_variables,
    BHEMeshData&& bhe_mesh_data)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data)),
      _bheMeshData(std::move(bhe_mesh_data)),
      _asm_mat_cache{_process_data._algebraic_BC_Setting._is_linear,
                     true /*use_monolithic_scheme*/}
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
    // BHE nodes with BHE type dependent number of variables.
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

        std::generate_n(std::back_inserter(all_mesh_subsets),
                        // Here the number of components equals to the
                        // number of unknowns on the BHE
                        number_of_unknowns,
                        [&ms = _mesh_subset_BHE_nodes.back()]()
                        { return *ms; });

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
        mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, element_to_bhe_map,
        mesh.isAxiallySymmetric(), _process_data);

    // Create BHE boundary conditions for each of the BHEs
    createBHEBoundaryConditionTopBottom(_bheMeshData.BHE_nodes);

    // Store BHE and soil elements to split the assembly and use the matrix
    // cache in the linear case only for soil elements
    if (_process_data._algebraic_BC_Setting._is_linear)
    {
        _bhes_element_ids = _bheMeshData.BHE_elements | ranges::views::join |
                            MeshLib::views::ids | ranges::to<std::vector>;

        // sort bhe elements if needed
        if (!std::is_sorted(_bhes_element_ids.begin(), _bhes_element_ids.end()))
        {
            std::sort(_bhes_element_ids.begin(), _bhes_element_ids.end());
        }

        _soil_element_ids = mesh.getElements() | MeshLib::views::ids |
                            ranges::to<std::vector>();

        // sort soil elements if needed
        if (!std::is_sorted(_soil_element_ids.begin(), _soil_element_ids.end()))
        {
            std::sort(_soil_element_ids.begin(), _soil_element_ids.end());
        }

        _soil_element_ids = ranges::views::set_difference(_soil_element_ids,
                                                          _bhes_element_ids) |
                            ranges::to<std::vector>();
    }

    if (_process_data._mass_lumping)
    {
        std::vector<std::size_t> const bhes_node_ids =
            _bheMeshData.BHE_nodes | ranges::views::join |
            ranges::views::transform([](auto const* const node)
                                     { return node->getID(); }) |
            ranges::to<std::vector>;

        // all connected soil elements and also the BHE elements.
        MeshLib::ElementSearch es{mesh};
        es.searchByNodeIDs(bhes_node_ids);

        assert(_process_data.mass_lumping_soil_elements.empty());
        _process_data.mass_lumping_soil_elements.resize(
            mesh.getNumberOfElements(), false);
        for (auto const id : es.getSearchedElementIDs())
        {
            _process_data.mass_lumping_soil_elements[id] = true;
        }
    }
}

void HeatTransportBHEProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble HeatTransportBHE process.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};

    if (_process_data._algebraic_BC_Setting._is_linear)
    {
        auto const& spec = this->getMatrixSpecifications(process_id);

        // use matrix cache for soil elements
        _asm_mat_cache.assemble(t, dt, x, x_prev, process_id, &M, &K, &b,
                                dof_table, _global_assembler, _local_assemblers,
                                _soil_element_ids);

        // reset the sparsity pattern for better performance in the BHE assembly
        MathLib::setMatrixSparsity(M, *spec.sparsity_pattern);
        MathLib::setMatrixSparsity(K, *spec.sparsity_pattern);

        // Call global assembler for each local BHE assembly item.
        GlobalExecutor::executeSelectedMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, _bhes_element_ids, dof_table, t, dt, x, x_prev,
            process_id, &M, &K, &b);
    }
    else
    {
        // Call global assembler for each local assembly item.
        GlobalExecutor::executeSelectedMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, getActiveElementIDs(), dof_table, t, dt, x,
            x_prev, process_id, &M, &K, &b);
    }

    // Algebraic BC procedure.
    if (_process_data._algebraic_BC_Setting._use_algebraic_bc)
    {
        algebraicBcConcreteProcess(t, dt, x, x_prev, process_id, M, K, b);
    }

    //_global_output(t, process_id, M, K, b);
}

void HeatTransportBHEProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HeatTransportBHE process.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, getActiveElementIDs(), dof_table, t, dt, x, x_prev,
        process_id, &b, &Jac);
}

void HeatTransportBHEProcess::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, int const process_id)
{
    DBUG("Compute heat flux for HeatTransportBHE process.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HeatTransportBHELocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, getActiveElementIDs(), dof_tables, t, dt, x, x_prev,
        process_id);
}

NumLib::IterationResult HeatTransportBHEProcess::postIterationConcreteProcess(
    GlobalVector const& x)
{
    // if the process use python boundary condition
    if (_process_data.py_bc_object == nullptr || !(_process_data._use_tespy))
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
        std::get<0>(_process_data.py_bc_object->dataframe_network),   // t
        std::get<1>(_process_data.py_bc_object->dataframe_network),   // T_in
        std::get<2>(_process_data.py_bc_object->dataframe_network));  // T_out
    if (!_process_data.py_bc_object->isOverriddenTespy())
    {
        DBUG("Method `tespySolver' not overridden in Python script.");
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

void HeatTransportBHEProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, const double t, const double dt,
    int const process_id)
{
    if (_process_data.py_bc_object == nullptr ||
        !_process_data._use_server_communication)
    {
        return;
    }

    auto& [time, Tin_value, Tout_value, Tout_nodes_ids, flowrate] =
        _process_data.py_bc_object->dataframe_network;

    // We found the problem that time != t, but it always equals the last
    // step. The following line is to correct this, although we do not use
    // it for server communication.
    time = t;

    auto const& solution = *x[process_id];

    // Iterate through each BHE
    const std::size_t n_bc_nodes = Tout_nodes_ids.size();
    for (std::size_t i = 0; i < n_bc_nodes; i++)
    {
        // read the T_out and store them in dataframe
        Tout_value[i] = solution[Tout_nodes_ids[i]];
    }

    // Transfer T_out to server_Communication and get back T_in and flowrate
    auto const server_communication_result =
        _process_data.py_bc_object->serverCommunicationPreTimestep(
            t, dt, Tin_value, Tout_value, flowrate);
    if (!_process_data.py_bc_object
             ->isOverriddenServerCommunicationPreTimestep())
    {
        DBUG("Method `serverCommunication' not overridden in Python script.");
    }

    auto const& [server_communication_Tin_value,
                 server_communication_flowrate] = server_communication_result;

    std::copy(begin(server_communication_Tin_value),
              end(server_communication_Tin_value),
              begin(Tin_value));
    std::copy(begin(server_communication_flowrate),
              end(server_communication_flowrate),
              begin(flowrate));
}

void HeatTransportBHEProcess::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& /*x_prev*/, const double t,
    const double dt, int const process_id)
{
    if (_process_data.py_bc_object == nullptr ||
        !_process_data._use_server_communication)
    {
        return;
    }

    auto& [time, Tin_value, Tout_value, Tout_nodes_ids, flowrate] =
        _process_data.py_bc_object->dataframe_network;

    // We found the problem that time != t, but it always equals the last
    // step. The following line is to correct this, although we do not use
    // it for server communication.
    time = t;

    auto const& solution = *x[process_id];

    // Iterate through each BHE
    const std::size_t n_bc_nodes = Tout_nodes_ids.size();
    for (std::size_t i = 0; i < n_bc_nodes; i++)
    {
        // read the T_out and store them in dataframe
        Tout_value[i] = solution[Tout_nodes_ids[i]];
    }

    // Transfer T_out to server_Communication
    _process_data.py_bc_object->serverCommunicationPostTimestep(
        t, dt, Tin_value, Tout_value, flowrate);
    if (!_process_data.py_bc_object
             ->isOverriddenServerCommunicationPostTimestep())
    {
        DBUG("Method `serverCommunication' not overridden in Python script.");
    }
}

void HeatTransportBHEProcess::algebraicBcConcreteProcess(
    [[maybe_unused]] const double t, double const /*dt*/,
    [[maybe_unused]] std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& /*xprev*/, int const /*process_id*/,
    [[maybe_unused]] GlobalMatrix& M, [[maybe_unused]] GlobalMatrix& K,
    [[maybe_unused]] GlobalVector& b)
{
#ifndef USE_PETSC
    auto M_normal = M.getRawMatrix();
    auto K_normal = K.getRawMatrix();
    auto n_original_rows = K_normal.rows();
    auto const n_BHE_bottom_pairs = _vec_bottom_BHE_node_indices.size();
    auto const n_BHE_top_pairs = _vec_top_BHE_node_indices.size();

    // apply weighting factor based on the max value from column wise inner
    // product and scale it with user defined value
    const double w_val =
        _process_data._algebraic_BC_Setting._weighting_factor *
        (Eigen::RowVectorXd::Ones(K_normal.rows()) * K_normal.cwiseAbs())
            .maxCoeff();

    M_normal.conservativeResize(
        M_normal.rows() + n_BHE_bottom_pairs + n_BHE_top_pairs,
        M_normal.cols());
    K_normal.conservativeResize(
        K_normal.rows() + n_BHE_bottom_pairs + n_BHE_top_pairs,
        K_normal.cols());

    for (std::size_t i = 0; i < n_BHE_bottom_pairs; i++)
    {
        Eigen::SparseVector<double> M_Plus(M_normal.cols());
        M_Plus.setZero();
        M_normal.row(n_original_rows + i) = M_Plus;

        Eigen::SparseVector<double> K_Plus(K_normal.cols());
        K_Plus.setZero();

        auto const [bhe_idx, first_BHE_bottom_index, second_BHE_bottom_index] =
            _vec_bottom_BHE_node_indices[i];

        K_Plus.insert(first_BHE_bottom_index) = w_val;
        K_Plus.insert(second_BHE_bottom_index) = -w_val;

        K_normal.row(n_original_rows + i) = K_Plus;
    }

    auto b_normal = b.getRawVector();
    Eigen::SparseVector<double> b_Plus(b_normal.rows() + n_BHE_bottom_pairs +
                                       n_BHE_top_pairs);
    b_Plus.setZero();

    // Copy values from the original column vector to the modified one
    for (int i = 0; i < b_normal.innerSize(); ++i)
    {
        b_Plus.insert(i) = b_normal.coeff(i);
    }

    for (std::size_t i = 0; i < n_BHE_top_pairs; i++)
    {
        Eigen::SparseVector<double> M_Plus(M_normal.cols());
        M_Plus.setZero();
        M_normal.row(n_original_rows + n_BHE_bottom_pairs + i) = M_Plus;

        Eigen::SparseVector<double> K_Plus(K_normal.cols());
        K_Plus.setZero();

        auto const [bhe_idx, first_BHE_top_index, second_BHE_top_index] =
            _vec_top_BHE_node_indices[i];

        auto first_BHE_top_index_pair = first_BHE_top_index;
        auto second_BHE_top_index_pair = second_BHE_top_index;

        K_Plus.insert(first_BHE_top_index_pair) =
            w_val;  // for power BC, the inflow node must be positive
        K_Plus.insert(second_BHE_top_index_pair) =
            -w_val;  // for power BC, the outflow node must be negative

        K_normal.row(n_original_rows + n_BHE_bottom_pairs + i) = K_Plus;

        // get the delta_T value here
        double const T_out = (*x[0])[second_BHE_top_index_pair];

        auto calculate_delta_T = [&](auto& bhe)
        {
            auto const T_in = bhe.updateFlowRateAndTemperature(T_out, t);
            return T_in - T_out;
        };
        auto delta_T = std::visit(calculate_delta_T,
                                  _process_data._vec_BHE_property[bhe_idx]);

        b_Plus.insert(n_original_rows + n_BHE_bottom_pairs + i) =
            delta_T * w_val;
    }

    M.getRawMatrix() = M_normal;
    K.getRawMatrix() = K_normal;
    b.getRawVector() = b_Plus;
#else
    OGS_FATAL(
        "The Algebraic Boundary Condition is not implemented for use with "
        "PETsc Library! Simulation will be terminated.");
#endif
}

void HeatTransportBHEProcess::createBHEBoundaryConditionTopBottom(
    std::vector<std::vector<MeshLib::Node*>> const& all_bhe_nodes)
{
    const int process_id = 0;
    auto& bcs = _boundary_conditions[process_id];

    std::size_t const n_BHEs = _process_data._vec_BHE_property.size();

    // for each BHE
    for (std::size_t bhe_i = 0; bhe_i < n_BHEs; bhe_i++)
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
            auto const& connected_elements =
                _mesh.getElementsConnectedToNode(*bhe_node);
            const std::size_t n_line_elements = std::count_if(
                connected_elements.begin(), connected_elements.end(),
                [](MeshLib::Element const* elem)
                { return (elem->getDimension() == 1); });

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
        if ((*bhe_boundary_nodes[0])[2] == (*bhe_boundary_nodes[1])[2])
        {
            INFO(
                "For 1P-type BHE, the BHE inflow and outflow "
                "nodes are identified according to their mesh node id in "
                "ascending order");
        }
        // for 1U, 2U, CXC, CXA type BHE
        else
        {
            // swap the boundary nodes if the z coordinate of the
            // first node is lower than it on the second node
            if ((*bhe_boundary_nodes[0])[2] < (*bhe_boundary_nodes[1])[2])
            {
                std::swap(bhe_boundary_nodes[0], bhe_boundary_nodes[1]);
            }
        }

        auto get_global_index =
            [&](std::size_t const node_id, int const component)
        {
            return _local_to_global_index_map->getGlobalIndex(
                {_mesh.getID(), MeshLib::MeshItemType::Node, node_id},
                variable_id, component);
        };

        auto get_global_bhe_bc_indices =
            [&](std::array<
                std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
                    nodes_and_components)
        {
            return std::make_pair(
                get_global_index(nodes_and_components[0].first,
                                 nodes_and_components[0].second),
                get_global_index(nodes_and_components[1].first,
                                 nodes_and_components[1].second));
        };

        auto get_global_bhe_bc_indices_with_bhe_idx =
            [&](std::size_t bhe_idx,
                std::array<
                    std::pair<std::size_t /*node_id*/, int /*component*/>, 2>
                    nodes_and_components)
        {
            return std::make_tuple(
                bhe_idx,
                get_global_index(nodes_and_components[0].first,
                                 nodes_and_components[0].second),
                get_global_index(nodes_and_components[1].first,
                                 nodes_and_components[1].second));
        };

        auto createBCs =
            [&, bc_top_node_id = bhe_boundary_nodes[0]->getID(),
             bc_bottom_node_id = bhe_boundary_nodes[1]->getID()](auto& bhe)
        {
            for (auto const& in_out_component_id :
                 bhe.inflow_outflow_bc_component_ids)
            {
                if (bhe.use_python_bcs ||
                    this->_process_data._use_server_communication)
                // call BHEPythonBoundarycondition
                {
                    if (this->_process_data
                            .py_bc_object)  // the bc object exist
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
                }
                else
                {
                    if (this->_process_data._algebraic_BC_Setting
                            ._use_algebraic_bc &&
                        bhe.isPowerBC())
                    {
                        // for algebraic_bc method, record the pair of indices
                        // in a separate vector
                        _vec_top_BHE_node_indices.push_back(
                            get_global_bhe_bc_indices_with_bhe_idx(
                                bhe_i,
                                {{{bc_top_node_id, in_out_component_id.first},
                                  {bc_top_node_id,
                                   in_out_component_id.second}}}));
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
                                    return bhe.updateFlowRateAndTemperature(T,
                                                                            t);
                                }));
                    }
                }

                auto const bottom_nodes_and_components =
                    bhe.getBHEBottomDirichletBCNodesAndComponents(
                        bc_bottom_node_id,
                        in_out_component_id.first,
                        in_out_component_id.second);

                if (bottom_nodes_and_components &&
                    !this->_process_data._algebraic_BC_Setting
                         ._use_algebraic_bc)
                {
                    // Bottom, outflow, all cases | not needed for algebraic_bc
                    // method
                    bcs.addBoundaryCondition(
                        createBHEBottomDirichletBoundaryCondition(
                            get_global_bhe_bc_indices(
                                {{{bc_bottom_node_id,
                                   in_out_component_id.first},
                                  {bc_bottom_node_id,
                                   in_out_component_id.second}}})));
                }
                else if (bottom_nodes_and_components &&
                         this->_process_data._algebraic_BC_Setting
                             ._use_algebraic_bc)
                {
                    // for algebraic_bc method, record the pair of indices in a
                    // separate vector
                    _vec_bottom_BHE_node_indices.push_back(
                        get_global_bhe_bc_indices_with_bhe_idx(
                            bhe_i,
                            {{{bc_bottom_node_id, in_out_component_id.first},
                              {bc_bottom_node_id,
                               in_out_component_id.second}}}));
                }
            }
        };
        visit(createBCs, _process_data._vec_BHE_property[bhe_i]);
    }
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
