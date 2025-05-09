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
template <int GlobalDim>
HydroMechanicsProcess<GlobalDim>::HydroMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HydroMechanicsProcessData<GlobalDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    INFO("[LIE/HM] looking for fracture elements in the given mesh");
    std::vector<int> vec_fracture_mat_IDs;
    std::vector<std::vector<MeshLib::Element*>> vec_vec_fracture_elements;
    std::vector<std::vector<MeshLib::Element*>>
        vec_vec_fracture_matrix_elements;
    std::vector<std::vector<MeshLib::Node*>> vec_vec_fracture_nodes;
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_branch_nodeID_matIDs;
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_junction_nodeID_matIDs;
    getFractureMatrixDataInMesh(
        mesh, _vec_matrix_elements, vec_fracture_mat_IDs,
        vec_vec_fracture_elements, vec_vec_fracture_matrix_elements,
        vec_vec_fracture_nodes, vec_branch_nodeID_matIDs,
        vec_junction_nodeID_matIDs);
    _vec_fracture_elements.insert(_vec_fracture_elements.begin(),
                                  vec_vec_fracture_elements[0].begin(),
                                  vec_vec_fracture_elements[0].end());
    _vec_fracture_matrix_elements.insert(
        _vec_fracture_matrix_elements.begin(),
        vec_vec_fracture_matrix_elements[0].begin(),
        vec_vec_fracture_matrix_elements[0].end());
    _vec_fracture_nodes.insert(_vec_fracture_nodes.begin(),
                               vec_vec_fracture_nodes[0].begin(),
                               vec_vec_fracture_nodes[0].end());

    if (!_vec_fracture_elements.empty())
    {
        // set fracture property assuming a fracture forms a straight line
        setFractureProperty(GlobalDim,
                            *_vec_fracture_elements[0],
                            *_process_data.fracture_property);
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
            if (std::find(vec_fracture_mat_IDs.begin(),
                          vec_fracture_mat_IDs.end(),
                          matID) == vec_fracture_mat_IDs.end())
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
}

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::constructDofTable()
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
    if (!_vec_fracture_nodes.empty())
    {
        // u jump
        _mesh_subset_fracture_nodes =
            std::make_unique<MeshLib::MeshSubset>(_mesh, _vec_fracture_nodes);
    }

    // Collect the mesh subsets in a vector.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::vector<int> vec_n_components;
    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    // pressure
    vec_n_components.push_back(1);
    all_mesh_subsets.emplace_back(*_mesh_subset_nodes_p);
    if (!_process_data.deactivate_matrix_in_flow)
    {
        vec_var_elements.push_back(&_mesh.getElements());
    }
    else
    {
        // TODO set elements including active nodes for pressure.
        // cannot use ElementStatus
        vec_var_elements.push_back(&_vec_fracture_matrix_elements);
    }
    // regular displacement
    vec_n_components.push_back(GlobalDim);
    std::generate_n(std::back_inserter(all_mesh_subsets), GlobalDim,
                    [&]() { return *_mesh_subset_matrix_nodes; });
    vec_var_elements.push_back(&_vec_matrix_elements);
    if (!_vec_fracture_nodes.empty())
    {
        // displacement jump
        vec_n_components.push_back(GlobalDim);
        std::generate_n(std::back_inserter(all_mesh_subsets), GlobalDim,
                        [&]() { return *_mesh_subset_fracture_nodes; });
        vec_var_elements.push_back(&_vec_fracture_matrix_elements);
    }

    INFO("[LIE/HM] creating a DoF table");
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    DBUG("created {:d} DoF", _local_to_global_index_map->size());
}

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    assert(mesh.getDimension() == GlobalDim);
    INFO("[LIE/HM] creating local assemblers");
    ProcessLib::LIE::HydroMechanics::createLocalAssemblers<
        GlobalDim, HydroMechanicsLocalAssemblerMatrix,
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

    add_secondary_variable(
        "sigma",
        MathLib::KelvinVector::KelvinVectorType<GlobalDim>::RowsAtCompileTime,
        &LocalAssemblerInterface::getIntPtSigma);

    add_secondary_variable(
        "epsilon",
        MathLib::KelvinVector::KelvinVectorType<GlobalDim>::RowsAtCompileTime,
        &LocalAssemblerInterface::getIntPtEpsilon);

    add_secondary_variable("velocity", GlobalDim,
                           &LocalAssemblerInterface::getIntPtDarcyVelocity);

    add_secondary_variable("fracture_velocity", GlobalDim,
                           &LocalAssemblerInterface::getIntPtFractureVelocity);

    add_secondary_variable("fracture_stress", GlobalDim,
                           &LocalAssemblerInterface::getIntPtFractureStress);

    add_secondary_variable("fracture_aperture", 1,
                           &LocalAssemblerInterface::getIntPtFractureAperture);

    add_secondary_variable(
        "fracture_permeability", 1,
        &LocalAssemblerInterface::getIntPtFracturePermeability);

    _process_data.element_stresses = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "sigma_avg",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::KelvinVectorType<GlobalDim>::RowsAtCompileTime);

    _process_data.element_velocities = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity_avg",
        MeshLib::MeshItemType::Cell, GlobalDim);

    if (!_vec_fracture_elements.empty())
    {
        auto mesh_prop_levelset = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "levelset1",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_levelset->resize(mesh.getNumberOfElements());
        for (MeshLib::Element const* e : _mesh.getElements())
        {
            if (e->getDimension() < GlobalDim)
            {
                continue;
            }

            std::vector<FractureProperty*> fracture_props(
                {_process_data.fracture_property.get()});
            std::vector<JunctionProperty*> junction_props;
            std::unordered_map<int, int> fracID_to_local({{0, 0}});
            std::vector<double> levelsets = uGlobalEnrichments(
                fracture_props, junction_props, fracID_to_local,
                Eigen::Vector3d(MeshLib::getCenterOfGravity(*e).data()));
            (*mesh_prop_levelset)[e->getID()] = levelsets[0];
        }

        _process_data.element_local_jumps =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "local_jump_w_avg",
                MeshLib::MeshItemType::Cell, GlobalDim);

        _process_data.element_fracture_stresses =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "fracture_stress_avg",
                MeshLib::MeshItemType::Cell, GlobalDim);

        _process_data.element_fracture_velocities =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "fracture_velocity_avg",
                MeshLib::MeshItemType::Cell, GlobalDim);

        auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "fracture_aperture_avg",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_b->resize(mesh.getNumberOfElements());

        auto const mesh_prop_matid = materialIDs(mesh);
        if (!mesh_prop_matid)
        {
            OGS_FATAL("Could not access MaterialIDs property from mesh.");
        }
        auto const& frac = _process_data.fracture_property;
        for (MeshLib::Element const* e : _mesh.getElements())
        {
            if (e->getDimension() == GlobalDim)
            {
                continue;
            }
            if ((*mesh_prop_matid)[e->getID()] != frac->mat_id)
            {
                continue;
            }
            // Mean value for the element. This allows usage of node based
            // properties for aperture.
            (*mesh_prop_b)[e->getID()] =
                frac->aperture0
                    .getNodalValuesOnElement(*e, /*time independent*/ 0)
                    .mean();
        }
        _process_data.mesh_prop_b = mesh_prop_b;

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
                MeshLib::MeshItemType::Node, GlobalDim);
        assert(_process_data.mesh_prop_nodal_forces->size() ==
               GlobalDim * mesh.getNumberOfNodes());

        _process_data.mesh_prop_nodal_forces_jump =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "NodalForcesJump",
                MeshLib::MeshItemType::Node, GlobalDim);
        assert(_process_data.mesh_prop_nodal_forces_jump->size() ==
               GlobalDim * mesh.getNumberOfNodes());

        _process_data.mesh_prop_hydraulic_flow =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "MassFlowRate",
                MeshLib::MeshItemType::Node, 1);
        assert(_process_data.mesh_prop_hydraulic_flow->size() ==
               mesh.getNumberOfNodes());
    }
}

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::postTimestepConcreteProcess(
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

template <int GlobalDim>
bool HydroMechanicsProcess<GlobalDim>::isLinear() const
{
    return false;
}

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::assembleConcreteProcess(
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

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::assembleWithJacobianConcreteProcess(
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

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::preTimestepConcreteProcess(
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
