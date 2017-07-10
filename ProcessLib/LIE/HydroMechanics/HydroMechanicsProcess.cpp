/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcess-fwd.h"

#include <algorithm>
#include <cassert>
#include <vector>

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/ElementStatus.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "ProcessLib/LIE/BoundaryCondition/BoundaryConditionBuilder.h"
#include "ProcessLib/LIE/Common/LevelSetFunction.h"
#include "ProcessLib/LIE/Common/MeshUtils.h"
#include "ProcessLib/LIE/Common/Utils.h"
#include "ProcessLib/LIE/HydroMechanics/LocalAssembler/HydroMechanicsLocalAssemblerFracture.h"
#include "ProcessLib/LIE/HydroMechanics/LocalAssembler/HydroMechanicsLocalAssemblerMatrix.h"
#include "ProcessLib/LIE/HydroMechanics/LocalAssembler/HydroMechanicsLocalAssemblerMatrixNearFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <unsigned GlobalDim>
HydroMechanicsProcess<GlobalDim>::HydroMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    HydroMechanicsProcessData<GlobalDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    INFO("[LIE/HM] looking for fracture elements in the given mesh");
    std::vector<int> vec_fracture_mat_IDs;
    std::vector<std::vector<MeshLib::Element*>> vec_vec_fracture_elements;
    std::vector<std::vector<MeshLib::Element*>>
        vec_vec_fracture_matrix_elements;
    std::vector<std::vector<MeshLib::Node*>> vec_vec_fracture_nodes;
    getFractureMatrixDataInMesh(mesh,
                                _vec_matrix_elements,
                                vec_fracture_mat_IDs,
                                vec_vec_fracture_elements,
                                vec_vec_fracture_matrix_elements,
                                vec_vec_fracture_nodes);
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
                            *_process_data.fracture_property.get());
    }

    // need to use a custom Neumann BC assembler for displacement jumps
    for (ProcessVariable& pv : getProcessVariables())
    {
        if (pv.getName().find("displacement_jump") == std::string::npos)
            continue;
        pv.setBoundaryConditionBuilder(
            std::make_unique<BoundaryConditionBuilder>(
                *_process_data.fracture_property.get()));
    }

    if (!_process_data.deactivate_matrix_in_flow)
    {
        _process_data.p_element_status =
            std::make_unique<MeshLib::ElementStatus>(&mesh);
    }
    else
    {
        auto range = MeshLib::MeshInformation::getValueBounds<int>(mesh, "MaterialIDs");
        std::vector<int> vec_p_inactive_matIDs;
        for (int matID = range.first; matID <= range.second; matID++)
        {
            if (std::find(vec_fracture_mat_IDs.begin(),
                          vec_fracture_mat_IDs.end(),
                          matID) == vec_fracture_mat_IDs.end())
                vec_p_inactive_matIDs.push_back(matID);
        }
        _process_data.p_element_status =
            std::make_unique<MeshLib::ElementStatus>(&mesh,
                                                     vec_p_inactive_matIDs);

        ProcessVariable const& pv_p = getProcessVariables()[0];
        _process_data.p0 = &pv_p.getInitialCondition();
    }
}


template <unsigned GlobalDim>
void HydroMechanicsProcess<GlobalDim>::constructDofTable()
{
    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // for extrapolation
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());
    // pressure
    _mesh_nodes_p = MeshLib::getBaseNodes(
        _process_data.p_element_status->getActiveElements());
    _mesh_subset_nodes_p =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh_nodes_p);
    // regular u
    _mesh_subset_matrix_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());
    if (!_vec_fracture_nodes.empty())
    {
        // u jump
        _mesh_subset_fracture_nodes =
            std::make_unique<MeshLib::MeshSubset>(_mesh, &_vec_fracture_nodes);
    }

    // Collect the mesh subsets in a vector.
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
    std::vector<unsigned> vec_n_components;
    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    // pressure
    vec_n_components.push_back(1);
    all_mesh_subsets.emplace_back(_mesh_subset_nodes_p.get());
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
    std::generate_n(std::back_inserter(all_mesh_subsets), GlobalDim, [&]() {
        return MeshLib::MeshSubsets{_mesh_subset_matrix_nodes.get()};
    });
    vec_var_elements.push_back(&_vec_matrix_elements);
    if (!_vec_fracture_nodes.empty())
    {
        // displacement jump
        vec_n_components.push_back(GlobalDim);
        std::generate_n(std::back_inserter(all_mesh_subsets), GlobalDim, [&]() {
            return MeshLib::MeshSubsets{_mesh_subset_fracture_nodes.get()};
        });
        vec_var_elements.push_back(&_vec_fracture_matrix_elements);
    }

    INFO("[LIE/HM] creating a DoF table");
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    DBUG("created %d DoF", _local_to_global_index_map->size());
}

template <unsigned GlobalDim>
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
        mesh.getElements(), dof_table,
        // use displacment process variable for shapefunction order
        getProcessVariables()[1].get().getShapeFunctionOrder(),
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

    auto mesh_prop_sigma_xx = MeshLib::getOrCreateMeshProperty<double>(
                                const_cast<MeshLib::Mesh&>(mesh),
                                "stress_xx", MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xx = mesh_prop_sigma_xx;

    auto mesh_prop_sigma_yy = MeshLib::getOrCreateMeshProperty<double>(
                                const_cast<MeshLib::Mesh&>(mesh),
                                "stress_yy", MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_yy = mesh_prop_sigma_yy;

    auto mesh_prop_sigma_xy = MeshLib::getOrCreateMeshProperty<double>(
                                const_cast<MeshLib::Mesh&>(mesh),
                                "stress_xy", MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xy = mesh_prop_sigma_xy;

    auto mesh_prop_epsilon_xx = MeshLib::getOrCreateMeshProperty<double>(
                                const_cast<MeshLib::Mesh&>(mesh),
                                "strain_xx", MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xx = mesh_prop_epsilon_xx;

    auto mesh_prop_epsilon_yy = MeshLib::getOrCreateMeshProperty<double>(
                                const_cast<MeshLib::Mesh&>(mesh),
                                "strain_yy", MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_yy = mesh_prop_epsilon_yy;

    auto mesh_prop_epsilon_xy = MeshLib::getOrCreateMeshProperty<double>(
                                const_cast<MeshLib::Mesh&>(mesh),
                                "strain_xy", MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xy = mesh_prop_epsilon_xy;

    auto mesh_prop_velocity = MeshLib::getOrCreateMeshProperty<double>(
                                const_cast<MeshLib::Mesh&>(mesh),
                                "velocity", MeshLib::MeshItemType::Cell, 3);
    mesh_prop_velocity->resize(mesh.getNumberOfElements() * 3);
    _process_data.mesh_prop_velocity = mesh_prop_velocity;

    if (!_vec_fracture_elements.empty())
    {
        auto mesh_prop_levelset = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "levelset1", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_levelset->resize(mesh.getNumberOfElements());
        for (MeshLib::Element const* e : _mesh.getElements())
        {
            if (e->getDimension() < GlobalDim)
                continue;

            double const levelsets =
                calculateLevelSetFunction(*_process_data.fracture_property,
                                          e->getCenterOfGravity().getCoords());
            (*mesh_prop_levelset)[e->getID()] = levelsets;
        }

        auto mesh_prop_w_n = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "w_n", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_w_n->resize(mesh.getNumberOfElements());
        auto mesh_prop_w_s = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "w_s", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_w_s->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_w_n = mesh_prop_w_n;
        _process_data.mesh_prop_w_s = mesh_prop_w_s;

        auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "aperture", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_b->resize(mesh.getNumberOfElements());
        auto mesh_prop_matid =
            mesh.getProperties().getPropertyVector<int>("MaterialIDs");
        auto frac = _process_data.fracture_property.get();
        for (MeshLib::Element const* e : _mesh.getElements())
        {
            if (e->getDimension() == GlobalDim)
                continue;
            if ((*mesh_prop_matid)[e->getID()] != frac->mat_id)
                continue;
            ProcessLib::SpatialPosition x;
            x.setElementID(e->getID());
            (*mesh_prop_b)[e->getID()] = (*frac->aperture0)(0, x)[0];
        }
        _process_data.mesh_prop_b = mesh_prop_b;


        auto mesh_prop_k_f = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "k_f", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_k_f->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_k_f = mesh_prop_k_f;

        auto mesh_prop_fracture_stress_shear = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "f_stress_s", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_stress_shear->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_shear =
            mesh_prop_fracture_stress_shear;

        auto mesh_prop_fracture_stress_normal = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "f_stress_n", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_stress_normal->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_normal =
            mesh_prop_fracture_stress_normal;

        auto mesh_prop_fracture_shear_failure = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "f_shear_failure", MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_shear_failure->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_shear_failure =
            mesh_prop_fracture_shear_failure;

        auto mesh_prop_nodal_w = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "nodal_w", MeshLib::MeshItemType::Node, GlobalDim);
        mesh_prop_nodal_w->resize(mesh.getNumberOfNodes() * GlobalDim);
        _process_data.mesh_prop_nodal_w = mesh_prop_nodal_w;

        auto mesh_prop_nodal_b = MeshLib::getOrCreateMeshProperty<double>(
                                    const_cast<MeshLib::Mesh&>(mesh),
                                    "nodal_aperture", MeshLib::MeshItemType::Node, 1);
        mesh_prop_nodal_b->resize(mesh.getNumberOfNodes());
        _process_data.mesh_prop_nodal_b = mesh_prop_nodal_b;
    }
}

template <unsigned GlobalDim>
void HydroMechanicsProcess<GlobalDim>::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x,
    StaggeredCouplingTerm const& coupled_term)
{
    DBUG("Compute the secondary variables for HydroMechanicsProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
        &HydroMechanicsLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, *_local_to_global_index_map, t, x, coupled_term);

    // Copy displacement jumps in a solution vector to mesh property
    // Remark: the copy is required because mesh properties for primary variables are set
    // during output and are not ready yet when this function is called.
    int g_variable_id = 0;
    int g_global_component_offset = 0;
    {
        int global_component_offset_next = 0;
        int global_component_offset = 0;
        for (int variable_id = 0;
             variable_id < static_cast<int>(this->getProcessVariables().size());
             ++variable_id)
        {
            ProcessVariable& pv = this->getProcessVariables()[variable_id];
            int const n_components = pv.getNumberOfComponents();
            global_component_offset = global_component_offset_next;
            global_component_offset_next += n_components;
            if (pv.getName() != "displacement_jump1")
                continue;

            g_variable_id = variable_id;
            g_global_component_offset = global_component_offset;
            break;
        }
    }

#ifdef USE_PETSC
    // TODO It is also possible directly to copy the data for single process
    // variable to a mesh property. It needs a vector of global indices and
    // some PETSc magic to do so.
    std::vector<double> x_copy(x.getLocalSize() + x.getGhostSize());
#else
    std::vector<double> x_copy(x.size());
#endif
    x.copyValues(x_copy);

    ProcessVariable& pv_g = this->getProcessVariables()[g_variable_id];
    auto& mesh_prop_g = pv_g.getOrCreateMeshProperty();
    auto const num_comp = pv_g.getNumberOfComponents();
    for (int component_id = 0; component_id < num_comp; ++component_id)
    {
        auto const& mesh_subsets =
            _local_to_global_index_map->getMeshSubsets(g_variable_id,
                                                       component_id);
        for (auto const& mesh_subset : mesh_subsets)
        {
            auto const mesh_id = mesh_subset->getMeshID();
            for (auto const* node : mesh_subset->getNodes())
            {
                MeshLib::Location const l(
                    mesh_id, MeshLib::MeshItemType::Node, node->getID());

                auto const global_component_id = g_global_component_offset + component_id;
                auto const index =
                        _local_to_global_index_map->getLocalIndex(
                            l, global_component_id, x.getRangeBegin(),
                            x.getRangeEnd());

                mesh_prop_g[node->getID() * num_comp + component_id] =
                        x_copy[index];
            }
        }
    }

    // compute nodal w and aperture
    auto const& R = _process_data.fracture_property->R;
    MeshLib::PropertyVector<double>& vec_w = *_process_data.mesh_prop_nodal_w;
    MeshLib::PropertyVector<double>& vec_b = *_process_data.mesh_prop_nodal_b;

    Eigen::VectorXd g(GlobalDim), w(GlobalDim);
    for (MeshLib::Node const* node : _vec_fracture_nodes)
    {
        auto const node_id = node->getID();
        g.setZero();
        for (unsigned k=0; k<GlobalDim; k++)
            g[k] = mesh_prop_g[node_id*GlobalDim + k];

        w.noalias() = R * g;
        for (unsigned k=0; k<GlobalDim; k++)
            vec_w[node_id*GlobalDim + k] = w[k];

        ProcessLib::SpatialPosition x;
        x.setNodeID(node_id);
        vec_b[node_id] = w[GlobalDim==2 ? 1 : 2] + (*_process_data.fracture_property->aperture0)(0,x)[0];
    }
}

// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class HydroMechanicsProcess<2>;

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
