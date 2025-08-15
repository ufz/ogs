/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HydroMechanicsProcess.h"

#include <range/v3/view/join.hpp>

#include "LocalAssembler/CreateLocalAssemblers.h"
#include "LocalAssembler/HydroMechanicsLocalAssemblerFracture.h"
#include "LocalAssembler/HydroMechanicsLocalAssemblerMatrix.h"
#include "LocalAssembler/HydroMechanicsLocalAssemblerMatrixNearFracture.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/ElementStatus.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "MeshToolsLib/MeshInformation.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/MeshElementParameter.h"
#include "ProcessLib/LIE/Common/BranchProperty.h"
#include "ProcessLib/LIE/Common/JunctionProperty.h"
#include "ProcessLib/LIE/Common/MeshUtils.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <int DisplacementDim>
HydroMechanicsProcess<DisplacementDim>::HydroMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HydroMechanicsProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    INFO("[LIE/HM] looking for fracture elements in the given mesh");
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_branch_nodeID_matIDs;
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_junction_nodeID_matIDs;
    getFractureMatrixDataInMesh(mesh, _vec_matrix_elements,
                                _vec_fracture_mat_IDs, _vec_fracture_elements,
                                _vec_fracture_matrix_elements,
                                _vec_fracture_nodes, vec_branch_nodeID_matIDs,
                                vec_junction_nodeID_matIDs);

    if (_vec_fracture_mat_IDs.size() !=
        _process_data.fracture_properties.size())
    {
        OGS_FATAL(
            "The number of the given fracture properties ({:d}) are not "
            "consistent with the number of fracture groups in a mesh ({:d}).",
            _process_data.fracture_properties.size(),
            _vec_fracture_mat_IDs.size());
    }

    // create a map from a material ID to a fracture ID
    auto max_frac_mat_id = std::max_element(_vec_fracture_mat_IDs.begin(),
                                            _vec_fracture_mat_IDs.end());
    _process_data.map_materialID_to_fractureID.resize(*max_frac_mat_id + 1);
    for (unsigned i = 0; i < _vec_fracture_mat_IDs.size(); i++)
    {
        _process_data.map_materialID_to_fractureID[_vec_fracture_mat_IDs[i]] =
            i;
    }

    // create a table of connected fracture IDs for each element
    _process_data.vec_ele_connected_fractureIDs.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < _vec_fracture_matrix_elements.size(); i++)
    {
        for (auto e : _vec_fracture_matrix_elements[i])
        {
            _process_data.vec_ele_connected_fractureIDs[e->getID()].push_back(
                i);
        }
    }

    // set fracture property
    for (auto& fracture_prop : _process_data.fracture_properties)
    {
        // based on the 1st element assuming a fracture forms a straight line
        setFractureProperty(
            DisplacementDim,
            *_vec_fracture_elements[fracture_prop.fracture_id][0],
            fracture_prop);
    }

    // set branches
    for (auto const& [vec_branch_nodeID, matID] : vec_branch_nodeID_matIDs)
    {
        auto master_matId = matID[0];
        auto slave_matId = matID[1];
        auto& master_frac =
            _process_data.fracture_properties
                [_process_data.map_materialID_to_fractureID[master_matId]];
        auto& slave_frac =
            _process_data.fracture_properties
                [_process_data.map_materialID_to_fractureID[slave_matId]];

        master_frac.branches_master.push_back(createBranchProperty(
            *mesh.getNode(vec_branch_nodeID), master_frac, slave_frac));

        slave_frac.branches_slave.push_back(createBranchProperty(
            *mesh.getNode(vec_branch_nodeID), master_frac, slave_frac));
    }

    // set junctions
    transform(cbegin(vec_junction_nodeID_matIDs),
              cend(vec_junction_nodeID_matIDs),
              back_inserter(_vec_junction_nodes),
              [&](auto& vec_junction_nodeID_matID)
              {
                  return const_cast<MeshLib::Node*>(
                      _mesh.getNode(vec_junction_nodeID_matID.first));
              });

    for (std::size_t i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto const& material_ids = vec_junction_nodeID_matIDs[i].second;
        assert(material_ids.size() == 2);
        std::array<int, 2> fracture_ids{
            {_process_data.map_materialID_to_fractureID[material_ids[0]],
             _process_data.map_materialID_to_fractureID[material_ids[1]]}};

        _process_data.junction_properties.emplace_back(
            i, *mesh.getNode(vec_junction_nodeID_matIDs[i].first),
            fracture_ids);
    }

    // create a table of connected junction IDs for each element
    _process_data.vec_ele_connected_junctionIDs.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto id :
             mesh.getElementsConnectedToNode(*node) | MeshLib::views::ids)
        {
            _process_data.vec_ele_connected_junctionIDs[id].push_back(i);
        }
    }

    // create a table of junction node and connected elements
    _vec_junction_fracture_matrix_elements.resize(
        vec_junction_nodeID_matIDs.size());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto e : mesh.getElementsConnectedToNode(*node))
        {
            _vec_junction_fracture_matrix_elements[i].push_back(
                const_cast<MeshLib::Element*>(e));
        }
    }

    //
    // If Neumann BCs for the displacement_jump variable are required they need
    // special treatment because of the levelset function. The implementation
    // exists in the version 6.1.0 (e54815cc07ee89c81f953a4955b1c788595dd725)
    // and was removed due to lack of applications.
    //

    if (!_process_data.deactivate_matrix_in_flow)
    {
        _process_data.p_element_status =
            std::make_unique<MeshLib::ElementStatus>(&mesh);
    }
    else
    {
        auto const range =
            MeshToolsLib::MeshInformation::getValueBounds(*materialIDs(mesh));
        if (!range)
        {
            OGS_FATAL(
                "Could not get minimum/maximum ranges values for the "
                "MaterialIDs property in the mesh '{:s}'.",
                mesh.getName());
        }

        std::vector<int> vec_p_inactive_matIDs;
        for (int matID = range->first; matID <= range->second; matID++)
        {
            if (std::find(_vec_fracture_mat_IDs.begin(),
                          _vec_fracture_mat_IDs.end(),
                          matID) == _vec_fracture_mat_IDs.end())
            {
                vec_p_inactive_matIDs.push_back(matID);
            }
        }
        _process_data.p_element_status =
            std::make_unique<MeshLib::ElementStatus>(&mesh,
                                                     vec_p_inactive_matIDs);

        const int monolithic_process_id = 0;
        ProcessVariable const& pv_p =
            getProcessVariables(monolithic_process_id)[0];
        _process_data.p0 = &pv_p.getInitialCondition();
    }
    MeshLib::PropertyVector<int> const* material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));
    _process_data.mesh_prop_materialIDs = material_ids;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::constructDofTable()
{
    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // for extrapolation
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    // pressure
    _mesh_nodes_p = MeshLib::getBaseNodes(
        _process_data.p_element_status->getActiveElements());
    _mesh_subset_nodes_p =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh_nodes_p);
    // regular u
    _mesh_subset_matrix_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    // u jump
    for (unsigned i = 0; i < _vec_fracture_nodes.size(); i++)
    {
        _mesh_subset_fracture_nodes.push_back(
            std::make_unique<MeshLib::MeshSubset const>(
                _mesh, _vec_fracture_nodes[i]));
    }
    // enrichment for junctions
    _mesh_subset_junction_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _vec_junction_nodes);

    // Collect the mesh subsets in a vector. (pressure, displacement,
    // displacement_jump_fracture, and displacement_jump_junction)
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    all_mesh_subsets.emplace_back(*_mesh_subset_nodes_p);
    std::generate_n(std::back_inserter(all_mesh_subsets), DisplacementDim,
                    [&]() { return *_mesh_subset_matrix_nodes; });
    for (auto const& ms : _mesh_subset_fracture_nodes)
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        DisplacementDim,
                        [&]() { return *ms; });
    }
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    DisplacementDim,
                    [&]() { return *_mesh_subset_junction_nodes; });

    // The corresponding number of components for (pressure, displacement,
    // displacement_jump_fracture, and displacement_jump_junction).
    std::vector<int> vec_n_components;
    vec_n_components.push_back(1);  // pressure
    vec_n_components.insert(
        vec_n_components.end(),
        1 + _vec_fracture_mat_IDs.size() + _vec_junction_nodes.size(),
        DisplacementDim);  // all displacements

    auto const all_fracture_matrix_elements = _vec_fracture_matrix_elements |
                                              ranges::views::join |
                                              ranges::to<std::vector>();
    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    if (!_process_data.deactivate_matrix_in_flow)
    {
        vec_var_elements.push_back(&_mesh.getElements());
    }
    else
    {
        // TODO set elements including active nodes for pressure.
        // cannot use ElementStatus
        vec_var_elements.push_back(&all_fracture_matrix_elements);
    }
    vec_var_elements.push_back(&_vec_matrix_elements);
    for (unsigned i = 0; i < _vec_fracture_matrix_elements.size(); i++)
    {
        vec_var_elements.push_back(&_vec_fracture_matrix_elements[i]);
    }
    for (unsigned i = 0; i < _vec_junction_fracture_matrix_elements.size(); i++)
    {
        vec_var_elements.push_back(&_vec_junction_fracture_matrix_elements[i]);
    }

    INFO("[LIE/HM] creating a DoF table");
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    DBUG("[LIE/HM] created {:d} DoF", _local_to_global_index_map->size());
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    assert(mesh.getDimension() == DisplacementDim);
    INFO("[LIE/HM] creating local assemblers");
    ProcessLib::LIE::HydroMechanics::createLocalAssemblers<
        DisplacementDim, HydroMechanicsLocalAssemblerMatrix,
        HydroMechanicsLocalAssemblerMatrixNearFracture,
        HydroMechanicsLocalAssemblerFracture>(
        mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    auto add_secondary_variable = [&](std::string const& name,
                                      int const num_components,
                                      auto get_ip_values_function)
    {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             _local_assemblers,
                             std::move(get_ip_values_function)));
    };

    add_secondary_variable("sigma",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtSigma);

    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtEpsilon);

    add_secondary_variable("velocity", DisplacementDim,
                           &LocalAssemblerInterface::getIntPtDarcyVelocity);

    add_secondary_variable("fracture_velocity", DisplacementDim,
                           &LocalAssemblerInterface::getIntPtFractureVelocity);

    add_secondary_variable("fracture_stress", DisplacementDim,
                           &LocalAssemblerInterface::getIntPtFractureStress);

    add_secondary_variable("fracture_aperture", 1,
                           &LocalAssemblerInterface::getIntPtFractureAperture);

    add_secondary_variable(
        "fracture_permeability", 1,
        &LocalAssemblerInterface::getIntPtFracturePermeability);

    _process_data.element_stresses = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "sigma_avg",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime);

    _process_data.element_velocities = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity_avg",
        MeshLib::MeshItemType::Cell, DisplacementDim);

    for (MeshLib::Element const* e : _mesh.getElements())
    {
        if (e->getDimension() < DisplacementDim)
        {
            continue;
        }

        Eigen::Vector3d const pt(getCenterOfGravity(*e).asEigenVector3d());
        std::vector<FractureProperty*> e_fracture_props;
        std::unordered_map<int, int> e_fracID_to_local;
        unsigned tmpi = 0;
        for (auto fid : _process_data.vec_ele_connected_fractureIDs[e->getID()])
        {
            e_fracture_props.push_back(&_process_data.fracture_properties[fid]);
            e_fracID_to_local.insert({fid, tmpi++});
        }
        std::vector<JunctionProperty*> e_junction_props;
        std::unordered_map<int, int> e_juncID_to_local;
        tmpi = 0;
        for (auto fid : _process_data.vec_ele_connected_junctionIDs[e->getID()])
        {
            e_junction_props.push_back(&_process_data.junction_properties[fid]);
            e_juncID_to_local.insert({fid, tmpi++});
        }
        std::vector<double> const levelsets(uGlobalEnrichments(
            e_fracture_props, e_junction_props, e_fracID_to_local, pt));

        for (unsigned i = 0; i < e_fracture_props.size(); i++)
        {
            auto mesh_prop_levelset = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh),
                "levelset" +
                    std::to_string(e_fracture_props[i]->fracture_id + 1),
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_levelset->resize(mesh.getNumberOfElements());
            (*mesh_prop_levelset)[e->getID()] = levelsets[i];
        }
        for (unsigned i = 0; i < e_junction_props.size(); i++)
        {
            auto mesh_prop_levelset = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh),
                "levelset" +
                    std::to_string(e_junction_props[i]->junction_id + 1 +
                                   _process_data.fracture_properties.size()),
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_levelset->resize(mesh.getNumberOfElements());
            (*mesh_prop_levelset)[e->getID()] =
                levelsets[i + e_fracture_props.size()];
        }
    }

    _process_data.element_local_jumps =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "local_jump_w_avg",
            MeshLib::MeshItemType::Cell, DisplacementDim);

    _process_data.element_fracture_stresses =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "fracture_stress_avg",
            MeshLib::MeshItemType::Cell, DisplacementDim);

    _process_data.element_fracture_velocities =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "fracture_velocity_avg",
            MeshLib::MeshItemType::Cell, DisplacementDim);

    auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "fracture_aperture_avg",
        MeshLib::MeshItemType::Cell, 1);

    mesh_prop_b->resize(mesh.getNumberOfElements());
    auto const& mesh_prop_matid = *_process_data.mesh_prop_materialIDs;
    for (auto const& fracture_prop : _process_data.fracture_properties)
    {
        for (MeshLib::Element const* e : _mesh.getElements())
        {
            if (e->getDimension() == DisplacementDim)
            {
                continue;
            }
            if (mesh_prop_matid[e->getID()] != fracture_prop.mat_id)
            {
                continue;
            }
            // Mean value for the element. This allows usage of node based
            // properties for aperture.
            (*mesh_prop_b)[e->getID()] =
                fracture_prop.aperture0
                    .getNodalValuesOnElement(*e, /*time independent*/ 0)
                    .mean();
        }
    }

    auto mesh_prop_k_f = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "fracture_permeability_avg",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_k_f->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_k_f = mesh_prop_k_f;

    auto mesh_prop_fracture_shear_failure =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "f_shear_failure",
            MeshLib::MeshItemType::Cell, 1);
    mesh_prop_fracture_shear_failure->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_fracture_shear_failure =
        mesh_prop_fracture_shear_failure;

    auto mesh_prop_nodal_p = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_nodal_p->resize(mesh.getNumberOfNodes());
    _process_data.mesh_prop_nodal_p = mesh_prop_nodal_p;

    _process_data.mesh_prop_nodal_forces =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "NodalForces",
            MeshLib::MeshItemType::Node, DisplacementDim);
    assert(_process_data.mesh_prop_nodal_forces->size() ==
           DisplacementDim * mesh.getNumberOfNodes());

    _process_data.mesh_prop_nodal_forces_jump =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "NodalForcesJump",
            MeshLib::MeshItemType::Node, DisplacementDim);
    assert(_process_data.mesh_prop_nodal_forces_jump->size() ==
           DisplacementDim * mesh.getNumberOfNodes());

    _process_data.mesh_prop_hydraulic_flow =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "MassFlowRate",
            MeshLib::MeshItemType::Node, 1);
    assert(_process_data.mesh_prop_hydraulic_flow->size() ==
           mesh.getNumberOfNodes());
    _process_data.mesh_prop_b = mesh_prop_b;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, const double t, double const dt,
    int const process_id)
{
    if (process_id == 0)
    {
        DBUG("PostTimestep HydroMechanicsProcess.");

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &HydroMechanicsLocalAssemblerInterface::postTimestep,
            _local_assemblers, getActiveElementIDs(), getDOFTables(x.size()), x,
            x_prev, t, dt, process_id);
    }

    DBUG("Compute the secondary variables for HydroMechanicsProcess.");

    const auto& dof_table = getDOFTable(process_id);

    // Copy displacement jumps in a solution vector to mesh property
    // Remark: the copy is required because mesh properties for primary
    // variables are set during output and are not ready yet when this function
    // is called.
    int g_variable_id = 0;
    {
        const int monolithic_process_id = 0;
        auto const& pvs = getProcessVariables(monolithic_process_id);
        auto const it =
            std::find_if(pvs.begin(), pvs.end(),
                         [](ProcessVariable const& pv)
                         { return pv.getName() == "displacement_jump1"; });
        if (it == pvs.end())
        {
            OGS_FATAL(
                "Didn't find expected 'displacement_jump1' process variable.");
        }
        g_variable_id = static_cast<int>(std::distance(pvs.begin(), it));
    }

    MathLib::LinAlg::setLocalAccessibleVector(*x[process_id]);

    const int monolithic_process_id = 0;
    ProcessVariable& pv_g =
        this->getProcessVariables(monolithic_process_id)[g_variable_id];
    auto const num_comp = pv_g.getNumberOfGlobalComponents();
    auto& mesh_prop_g = *MeshLib::getOrCreateMeshProperty<double>(
        _mesh, pv_g.getName(), MeshLib::MeshItemType::Node, num_comp);
    for (int component_id = 0; component_id < num_comp; ++component_id)
    {
        auto const& mesh_subset =
            dof_table.getMeshSubset(g_variable_id, component_id);
        for (auto const& l : MeshLib::views::meshLocations(
                 mesh_subset, MeshLib::MeshItemType::Node))
        {
            auto const global_index =
                dof_table.getGlobalIndex(l, g_variable_id, component_id);
            auto const node_id = l.item_id;
            mesh_prop_g[node_id * num_comp + component_id] =
                (*x[process_id])[global_index];
        }
    }
}

template <int DisplacementDim>
bool HydroMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble HydroMechanicsProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, dt, x, x_prev, process_id, &M, &K, &b);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HydroMechanicsProcess.");

    // Call global assembler for each local assembly item.
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, getActiveElementIDs(), dof_table, t, dt, x, x_prev,
        process_id, &b, &Jac);

    auto copyRhs = [&](int const variable_id, auto& output_vector)
    {
        transformVariableFromGlobalVector(b, variable_id,
                                          *_local_to_global_index_map,
                                          output_vector, std::negate<double>());
    };
    copyRhs(0, *_process_data.mesh_prop_hydraulic_flow);
    copyRhs(1, *_process_data.mesh_prop_nodal_forces);
    copyRhs(2, *_process_data.mesh_prop_nodal_forces_jump);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    int const process_id)
{
    DBUG("PreTimestep HydroMechanicsProcess.");

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HydroMechanicsLocalAssemblerInterface::preTimestep, _local_assemblers,
        getActiveElementIDs(), *_local_to_global_index_map, *x[process_id], t,
        dt);
}

// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class HydroMechanicsProcess<2>;
template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
