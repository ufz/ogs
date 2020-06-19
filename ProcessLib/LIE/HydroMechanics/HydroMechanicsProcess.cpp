/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HydroMechanicsProcess.h"

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/ElementStatus.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Properties.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "ParameterLib/MeshElementParameter.h"
#include "ProcessLib/LIE/Common/BranchProperty.h"
#include "ProcessLib/LIE/Common/JunctionProperty.h"
#include "ProcessLib/LIE/Common/MeshUtils.h"

#include "LocalAssembler/CreateLocalAssemblers.h"
#include "LocalAssembler/HydroMechanicsLocalAssemblerFracture.h"
#include "LocalAssembler/HydroMechanicsLocalAssemblerMatrix.h"
#include "LocalAssembler/HydroMechanicsLocalAssemblerMatrixNearFracture.h"

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
            MeshLib::MeshInformation::getValueBounds<int>(mesh, "MaterialIDs");
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
    const int monolithic_process_id = 0;
    ProcessLib::LIE::HydroMechanics::createLocalAssemblers<
        GlobalDim, HydroMechanicsLocalAssemblerMatrix,
        HydroMechanicsLocalAssemblerMatrixNearFracture,
        HydroMechanicsLocalAssemblerFracture>(
        mesh.getElements(), dof_table,
        // use displacment process variable for shapefunction order
        getProcessVariables(
            monolithic_process_id)[1].get().getShapeFunctionOrder(),
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);

    auto mesh_prop_sigma_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xx = mesh_prop_sigma_xx;

    auto mesh_prop_sigma_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_yy = mesh_prop_sigma_yy;

    auto mesh_prop_sigma_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_zz->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_zz = mesh_prop_sigma_zz;

    auto mesh_prop_sigma_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xy = mesh_prop_sigma_xy;

    if (GlobalDim == 3)
    {
        auto mesh_prop_sigma_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_xz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_stress_xz = mesh_prop_sigma_xz;

        auto mesh_prop_sigma_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_yz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_stress_yz = mesh_prop_sigma_yz;
    }

    auto mesh_prop_epsilon_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xx = mesh_prop_epsilon_xx;

    auto mesh_prop_epsilon_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_yy = mesh_prop_epsilon_yy;

    auto mesh_prop_epsilon_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_zz->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_zz = mesh_prop_epsilon_zz;

    auto mesh_prop_epsilon_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xy = mesh_prop_epsilon_xy;

    if (GlobalDim == 3)
    {
        auto mesh_prop_epsilon_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_xz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_strain_xz = mesh_prop_epsilon_xz;

        auto mesh_prop_epsilon_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_yz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_strain_yz = mesh_prop_epsilon_yz;
    }

    auto mesh_prop_velocity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity",
        MeshLib::MeshItemType::Cell, 3);
    mesh_prop_velocity->resize(mesh.getNumberOfElements() * 3);
    _process_data.mesh_prop_velocity = mesh_prop_velocity;

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
                Eigen::Vector3d(e->getCenterOfGravity().getCoords()));
            (*mesh_prop_levelset)[e->getID()] = levelsets[0];
        }

        auto mesh_prop_w_n = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "w_n",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_w_n->resize(mesh.getNumberOfElements());
        auto mesh_prop_w_s = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "w_s",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_w_s->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_w_n = mesh_prop_w_n;
        _process_data.mesh_prop_w_s = mesh_prop_w_s;

        auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "aperture",
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
            const_cast<MeshLib::Mesh&>(mesh), "k_f",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_k_f->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_k_f = mesh_prop_k_f;

        auto mesh_prop_fracture_stress_shear =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "f_stress_s",
                MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_stress_shear->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_shear =
            mesh_prop_fracture_stress_shear;

        auto mesh_prop_fracture_stress_normal =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "f_stress_n",
                MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_stress_normal->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_normal =
            mesh_prop_fracture_stress_normal;

        auto mesh_prop_fracture_shear_failure =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "f_shear_failure",
                MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_shear_failure->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_shear_failure =
            mesh_prop_fracture_shear_failure;

        auto mesh_prop_nodal_w = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "nodal_w",
            MeshLib::MeshItemType::Node, GlobalDim);
        mesh_prop_nodal_w->resize(mesh.getNumberOfNodes() * GlobalDim);
        _process_data.mesh_prop_nodal_w = mesh_prop_nodal_w;

        auto mesh_prop_nodal_b = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "nodal_aperture",
            MeshLib::MeshItemType::Node, 1);
        mesh_prop_nodal_b->resize(mesh.getNumberOfNodes());
        _process_data.mesh_prop_nodal_b = mesh_prop_nodal_b;

        if (GlobalDim == 3)
        {
            auto mesh_prop_w_s2 = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "w_s2",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_w_s2->resize(mesh.getNumberOfElements());
            _process_data.mesh_prop_w_s2 = mesh_prop_w_s2;

            auto mesh_prop_fracture_stress_shear2 =
                MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "f_stress_s2",
                    MeshLib::MeshItemType::Cell, 1);
            mesh_prop_fracture_stress_shear2->resize(
                mesh.getNumberOfElements());
            _process_data.mesh_prop_fracture_stress_shear2 =
                mesh_prop_fracture_stress_shear2;
        }

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
                const_cast<MeshLib::Mesh&>(mesh), "HydraulicFlow",
                MeshLib::MeshItemType::Node, 1);
        assert(_process_data.mesh_prop_hydraulic_flow->size() ==
               mesh.getNumberOfNodes());
    }
}

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, const double t, double const dt,
    int const process_id)
{
    DBUG("Compute the secondary variables for HydroMechanicsProcess.");
    const auto& dof_table = getDOFTable(process_id);

    {
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &HydroMechanicsLocalAssemblerInterface::postTimestep,
            _local_assemblers, pv.getActiveElementIDs(), dof_table,
            *x[process_id], t, dt);
    }

    // Copy displacement jumps in a solution vector to mesh property
    // Remark: the copy is required because mesh properties for primary
    // variables are set during output and are not ready yet when this function
    // is called.
    int g_variable_id = 0;
    {
        const int monolithic_process_id = 0;
        auto const& pvs = getProcessVariables(monolithic_process_id);
        auto const it =
            std::find_if(pvs.begin(), pvs.end(), [](ProcessVariable const& pv) {
                return pv.getName() == "displacement_jump1";
            });
        if (it == pvs.end())
        {
            OGS_FATAL(
                "Didn't find expected 'displacement_jump1' process "
                "variable.");
        }
        g_variable_id = static_cast<int>(std::distance(pvs.begin(), it));
    }

    MathLib::LinAlg::setLocalAccessibleVector(*x[process_id]);

    const int monolithic_process_id = 0;
    ProcessVariable& pv_g =
        this->getProcessVariables(monolithic_process_id)[g_variable_id];
    auto const num_comp = pv_g.getNumberOfComponents();
    auto& mesh_prop_g = *MeshLib::getOrCreateMeshProperty<double>(
        _mesh, pv_g.getName(), MeshLib::MeshItemType::Node, num_comp);
    for (int component_id = 0; component_id < num_comp; ++component_id)
    {
        auto const& mesh_subset = dof_table.getMeshSubset(
            g_variable_id, component_id);
        auto const mesh_id = mesh_subset.getMeshID();
        for (auto const* node : mesh_subset.getNodes())
        {
            MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                      node->getID());

            auto const global_index =
                dof_table.getGlobalIndex(l, g_variable_id, component_id);
            mesh_prop_g[node->getID() * num_comp + component_id] =
                (*x[process_id])[global_index];
        }
    }

    // compute nodal w and aperture
    auto const& R = _process_data.fracture_property->R;
    auto const& b0 = _process_data.fracture_property->aperture0;
    MeshLib::PropertyVector<double>& vec_w = *_process_data.mesh_prop_nodal_w;
    MeshLib::PropertyVector<double>& vec_b = *_process_data.mesh_prop_nodal_b;

    auto compute_nodal_aperture = [&](std::size_t const node_id,
                                      double const w_n) {
        // skip aperture computation for element-wise defined b0 because there
        // are jumps on the nodes between the element's values.
        if (dynamic_cast<ParameterLib::MeshElementParameter<double> const*>(
                &b0))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        ParameterLib::SpatialPosition x;
        x.setNodeID(node_id);
        return w_n + b0(/*time independent*/ 0, x)[0];
    };

    Eigen::VectorXd g(GlobalDim);
    Eigen::VectorXd w(GlobalDim);
    for (MeshLib::Node const* node : _vec_fracture_nodes)
    {
        auto const node_id = node->getID();
        g.setZero();
        for (int k = 0; k < GlobalDim; k++)
        {
            g[k] = mesh_prop_g[node_id * GlobalDim + k];
        }

        w.noalias() = R * g;
        for (int k = 0; k < GlobalDim; k++)
        {
            vec_w[node_id * GlobalDim + k] = w[k];
        }

        vec_b[node_id] = compute_nodal_aperture(node_id, w[GlobalDim - 1]);
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
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble HydroMechanicsProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, dt, x, xdot, process_id, M, K, b, _coupled_solutions);
}

template <int GlobalDim>
void HydroMechanicsProcess<GlobalDim>::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HydroMechanicsProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);

    auto copyRhs = [&](int const variable_id, auto& output_vector) {
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

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HydroMechanicsLocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map, *x[process_id],
        t, dt);
}

// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class HydroMechanicsProcess<2>;
template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
