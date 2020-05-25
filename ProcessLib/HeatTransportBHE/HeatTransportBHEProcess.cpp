/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
      process_data_(std::move(process_data)),
      bheMeshData_(getBHEDataInMesh(mesh))
{
    if (bheMeshData_.BHE_mat_IDs.size() !=
        process_data_.vec_BHE_property_.size())
    {
        OGS_FATAL(
            "The number of the given BHE properties ({:d}) are not consistent "
            "with the number of BHE groups in the mesh ({:d}).",
            process_data_.vec_BHE_property_.size(),
            bheMeshData_.BHE_mat_IDs.size());
    }

    auto material_ids = MeshLib::materialIDs(mesh);
    if (material_ids == nullptr)
    {
        OGS_FATAL("Not able to get material IDs! ");
    }

    process_data_.mesh_prop_materialIDs_ = material_ids;

    // create a map from a material ID to a BHE ID
    for (int i = 0; i < static_cast<int>(bheMeshData_.BHE_mat_IDs.size()); i++)
    {
        // fill in the map structure
        process_data_.map_materialID_to_BHE_ID_[bheMeshData_.BHE_mat_IDs[i]] =
            i;
    }
}

void HeatTransportBHEProcess::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    mesh_subset_all_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, mesh_.getNodes());

    //
    // Soil temperature variable defined on the whole mesh.
    //
    mesh_subset_soil_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, mesh_.getNodes());
    std::vector<MeshLib::MeshSubset> all_mesh_subsets{*mesh_subset_soil_nodes_};

    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    vec_var_elements.push_back(&(mesh_.getElements()));

    std::vector<int> vec_n_components{
        1};  // one component for the soil temperature variable.

    //
    // BHE nodes with BHE type dependend number of variables.
    //
    int const n_BHEs = process_data_.vec_BHE_property_.size();
    assert(n_BHEs == static_cast<int>(bheMeshData_.BHE_mat_IDs.size()));
    assert(n_BHEs == static_cast<int>(bheMeshData_.BHE_nodes.size()));
    assert(n_BHEs == static_cast<int>(bheMeshData_.BHE_elements.size()));

    // the BHE nodes need to be cherry-picked from the vector
    for (int i = 0; i < n_BHEs; i++)
    {
        auto const number_of_unknowns =
            visit([](auto const& bhe) { return bhe.number_of_unknowns; },
                  process_data_.vec_BHE_property_[i]);
        auto const& bhe_nodes = bheMeshData_.BHE_nodes[i];
        auto const& bhe_elements = bheMeshData_.BHE_elements[i];

        // All the BHE nodes have additional variables.
        mesh_subset_BHE_nodes_.push_back(
            std::make_unique<MeshLib::MeshSubset const>(mesh_, bhe_nodes));

        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            // Here the number of components equals to the
            // number of unknowns on the BHE
            number_of_unknowns,
            [& ms = mesh_subset_BHE_nodes_.back()]() { return *ms; });

        vec_n_components.push_back(number_of_unknowns);
        vec_var_elements.push_back(&bhe_elements);
    }

    local_to_global_index_map_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    // in case of debugging the dof table, activate the following line
    // std::cout << *local_to_global_index_map_ << "\n";
}

void HeatTransportBHEProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    // Quick access map to BHE's through element ids.
    std::unordered_map<std::size_t, BHE::BHETypes*> element_to_bhe_map;
    int const n_BHEs = process_data_.vec_BHE_property_.size();
    for (int i = 0; i < n_BHEs; i++)
    {
        auto const& bhe_elements = bheMeshData_.BHE_elements[i];
        for (auto const& e : bhe_elements)
        {
            element_to_bhe_map[e->getID()] =
                &process_data_.vec_BHE_property_[i];
        }
    }

    assert(mesh.getDimension() == 3);
    ProcessLib::HeatTransportBHE::createLocalAssemblers<
        HeatTransportBHELocalAssemblerSoil, HeatTransportBHELocalAssemblerBHE>(
        mesh.getElements(), dof_table, local_assemblers_, element_to_bhe_map,
        mesh.isAxiallySymmetric(), integration_order, process_data_);

    // Create BHE boundary conditions for each of the BHEs
    createBHEBoundaryConditionTopBottom(bheMeshData_.BHE_nodes);
}

void HeatTransportBHEProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble HeatTransportBHE process.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

void HeatTransportBHEProcess::assembleWithJacobianConcreteProcess(
    const double /*t*/, double const /*dt*/,
    std::vector<GlobalVector*> const& /*x*/, GlobalVector const& /*xdot*/,
    const double /*dxdot_dx*/, const double /*dx_dx*/, int const /*process_id*/,
    GlobalMatrix& /*M*/, GlobalMatrix& /*K*/, GlobalVector& /*b*/,
    GlobalMatrix& /*Jac*/)
{
    OGS_FATAL(
        "HeatTransportBHE: analytical Jacobian assembly is not implemented");
}

void HeatTransportBHEProcess::computeSecondaryVariableConcrete(
    double const t, double const dt, GlobalVector const& x,
    GlobalVector const& x_dot, int const process_id)
{
    DBUG("Compute heat flux for HeatTransportBHE process.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HeatTransportBHELocalAssemblerInterface::computeSecondaryVariable,
        local_assemblers_, pv.getActiveElementIDs(), getDOFTable(process_id), t,
        dt, x, x_dot, coupled_solutions_);
}

#ifdef OGS_USE_PYTHON
NumLib::IterationResult HeatTransportBHEProcess::postIterationConcreteProcess(
    GlobalVector const& x)
{
    // if the process use python boundary conditon
    if (process_data_.py_bc_object == nullptr)
        return NumLib::IterationResult::SUCCESS;

    // Here the task is to get current time flowrate and flow temperature from
    // TESPy
    auto const Tout_nodes_id =
        std::get<3>(process_data_.py_bc_object->dataframe_network);
    const std::size_t n_bc_nodes = Tout_nodes_id.size();

    // update flowrate in network if network exist a dynamic flowrate in time
    auto const cur_time =
        std::get<0>(process_data_.py_bc_object->dataframe_network);
    if (std::get<0>(process_data_.py_bc_object->tespyHydroSolver(cur_time)))
    {
        // calculate the current flowrate in each BHE from TESPy
        auto const cur_flowrate =
            std::get<1>(process_data_.py_bc_object->tespyHydroSolver(cur_time));
        for (std::size_t i = 0; i < n_bc_nodes; i++)
            std::get<4>(process_data_.py_bc_object->dataframe_network)[i] =
                cur_flowrate[i];
        if (!process_data_.py_bc_object->isOverriddenTespyHydro())
        {
            DBUG(
                "Method `tespyHydroSolver' not overridden in Python "
                "script.");
        }
    }

    // get the outflow temperature,
    // transfer it to TESPy and to get inflow temperature,
    // and determine whether it converges.
    for (std::size_t i = 0; i < n_bc_nodes; i++)
    {
        // read the T_out and store them in dataframe
        std::get<2>(process_data_.py_bc_object->dataframe_network)[i] =
            x[Tout_nodes_id[i]];
    }
    // Tout transfer to Python
    auto const tespy_result = process_data_.py_bc_object->tespyThermalSolver(
        std::get<0>(process_data_.py_bc_object->dataframe_network),
        std::get<1>(process_data_.py_bc_object->dataframe_network),
        std::get<2>(process_data_.py_bc_object->dataframe_network));
    if (!process_data_.py_bc_object->isOverriddenTespyThermal())
    {
        DBUG(
            "Method `tespyThermalSolver' not overridden in Python "
            "script.");
    }
    auto const cur_Tin = std::get<2>(tespy_result);

    // update the T_in
    for (std::size_t i = 0; i < n_bc_nodes; i++)
        std::get<1>(process_data_.py_bc_object->dataframe_network)[i] =
            cur_Tin[i];

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
    auto& bcs = boundary_conditions_[process_id];

    int const n_BHEs = static_cast<int>(process_data_.vec_BHE_property_.size());

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
        auto get_global_index =
            [&](std::size_t const node_id, int const component) {
                return local_to_global_index_map_->getGlobalIndex(
                    {mesh_.getID(), MeshLib::MeshItemType::Node, node_id},
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
                    if (process_data_.py_bc_object)  // the bc object exist
                    {
                        // apply the customized top, inflow BC.
                        bcs.addBoundaryCondition(
                            ProcessLib::createBHEInflowPythonBoundaryCondition(
                                get_global_bhe_bc_indices(
                                    bhe.getBHEInflowDirichletBCNodesAndComponents(
                                        bc_top_node_id, bc_bottom_node_id,
                                        in_out_component_id.first)),
                                bhe,
                                *(process_data_.py_bc_object)));
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
        visit(createBCs, process_data_.vec_BHE_property_[bhe_i]);
    }
}
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
