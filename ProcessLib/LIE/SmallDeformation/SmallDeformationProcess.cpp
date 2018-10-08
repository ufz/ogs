/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationProcess.h"

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "ProcessLib/LIE/Common/BranchProperty.h"
#include "ProcessLib/LIE/Common/MeshUtils.h"
#include "ProcessLib/LIE/SmallDeformation/LocalAssembler/CreateLocalAssemblers.h"
#include "ProcessLib/LIE/SmallDeformation/LocalAssembler/SmallDeformationLocalAssemblerFracture.h"
#include "ProcessLib/LIE/SmallDeformation/LocalAssembler/SmallDeformationLocalAssemblerMatrix.h"
#include "ProcessLib/LIE/SmallDeformation/LocalAssembler/SmallDeformationLocalAssemblerMatrixNearFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <int DisplacementDim>
SmallDeformationProcess<DisplacementDim>::SmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
        jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order,
              std::move(process_variables),
              std::move(secondary_variables),
              std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    std::vector<std::pair<std::size_t,std::vector<int>>> vec_branch_nodeID_matIDs;
    std::vector<std::pair<std::size_t,std::vector<int>>> vec_junction_nodeID_matIDs;
    getFractureMatrixDataInMesh(mesh,
                                _vec_matrix_elements,
                                _vec_fracture_mat_IDs,
                                _vec_fracture_elements,
                                _vec_fracture_matrix_elements,
                                _vec_fracture_nodes,
                                vec_branch_nodeID_matIDs, vec_junction_nodeID_matIDs);

    if (_vec_fracture_mat_IDs.size() !=
        _process_data._vec_fracture_property.size())
    {
        OGS_FATAL(
            "The number of the given fracture properties (%d) are not "
            "consistent"
            " with the number of fracture groups in a mesh (%d).",
            _process_data._vec_fracture_property.size(),
            _vec_fracture_mat_IDs.size());
    }

    // create a map from a material ID to a fracture ID
    auto max_frac_mat_id = std::max_element(_vec_fracture_mat_IDs.begin(),
                                            _vec_fracture_mat_IDs.end());
    _process_data._map_materialID_to_fractureID.resize(*max_frac_mat_id + 1);
    for (unsigned i = 0; i < _vec_fracture_mat_IDs.size(); i++)
    {
        _process_data._map_materialID_to_fractureID[_vec_fracture_mat_IDs[i]] =
            i;
    }

    // create a table of connected fracture IDs for each element
    _process_data._vec_ele_connected_fractureIDs.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < _vec_fracture_matrix_elements.size(); i++)
    {
        for (auto e : _vec_fracture_matrix_elements[i])
        {
            _process_data._vec_ele_connected_fractureIDs[e->getID()].push_back(
                i);
        }
    }

    // set fracture property
    for (auto& fracture_prop : _process_data._vec_fracture_property)
    {
        // based on the 1st element assuming a fracture forms a straight line
        setFractureProperty(
            DisplacementDim,
            *_vec_fracture_elements[fracture_prop->fracture_id][0],
            *fracture_prop.get());
    }

    // set branches
    for (std::size_t i=0; i<vec_branch_nodeID_matIDs.size(); i++)
    {
        auto master_matId = vec_branch_nodeID_matIDs[i].second[0];
        auto slave_matId = vec_branch_nodeID_matIDs[i].second[1];
        auto& master_frac =
            *_process_data._vec_fracture_property
                 [_process_data._map_materialID_to_fractureID[master_matId]];
        auto& slave_frac =
            *_process_data._vec_fracture_property
                 [_process_data._map_materialID_to_fractureID[slave_matId]];

        BranchProperty* branch = new BranchProperty();
        setBranchProperty(*mesh.getNode(vec_branch_nodeID_matIDs[i].first),
                          master_frac, slave_frac, *branch);

        master_frac.branches_master.emplace_back(branch);

        BranchProperty* branch2 = new BranchProperty();
        setBranchProperty(*mesh.getNode(vec_branch_nodeID_matIDs[i].first),
                          master_frac, slave_frac, *branch2);
        slave_frac.branches_slave.emplace_back(branch2);
    }

    // set junctions
    for (std::size_t i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        _vec_junction_nodes.push_back(
            const_cast<MeshLib::Node*>(_mesh.getNode(vec_junction_nodeID_matIDs[i].first)));
    }
    for (std::size_t i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        JunctionProperty* junction = new JunctionProperty();
        std::vector<int> fracIDs;
        for (auto matid : vec_junction_nodeID_matIDs[i].second)
            fracIDs.push_back(_process_data._map_materialID_to_fractureID[matid]);
        setJunctionProperty(i, *mesh.getNode(vec_junction_nodeID_matIDs[i].first), fracIDs,
                            *junction);

        _process_data._vec_junction_property.emplace_back(junction);
    }

    // create a table of connected junction IDs for each element
    _process_data._vec_ele_connected_junctionIDs.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto e : node->getElements())
        {
            _process_data._vec_ele_connected_junctionIDs[e->getID()].push_back(
                i);
        }
    }

    // create a table of junction node and connected elements
    _vec_junction_fracture_matrix_elements.resize(
        vec_junction_nodeID_matIDs.size());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto e : node->getElements())
        {
            _vec_junction_fracture_matrix_elements[i].push_back(e);
        }
    }

    //
    // If Neumann BCs for the displacement_jump variable are required they need
    // special treatment because of the levelset function. The implementation
    // exists in the version 6.1.0 (e54815cc07ee89c81f953a4955b1c788595dd725)
    // and was removed due to lack of applications.
    //

    MeshLib::PropertyVector<int> const* material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));
    _process_data._mesh_prop_materialIDs = material_ids;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::constructDofTable()
{
    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // for extrapolation
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
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

    // Collect the mesh subsets in a vector.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::generate_n(std::back_inserter(all_mesh_subsets), DisplacementDim,
                    [&]() { return *_mesh_subset_matrix_nodes; });
    for (auto& ms : _mesh_subset_fracture_nodes)
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        DisplacementDim,
                        [&]() { return *ms; });
    }
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    DisplacementDim,
                    [&]() { return *_mesh_subset_junction_nodes; });

    std::vector<int> const vec_n_components(1 + _vec_fracture_mat_IDs.size() + _vec_junction_nodes.size(),
                                            DisplacementDim);

    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    vec_var_elements.push_back(&_vec_matrix_elements);
    for (unsigned i = 0; i < _vec_fracture_matrix_elements.size(); i++)
    {
        vec_var_elements.push_back(&_vec_fracture_matrix_elements[i]);
    }
    for (unsigned i = 0; i < _vec_junction_fracture_matrix_elements.size(); i++)
    {
        vec_var_elements.push_back(&_vec_junction_fracture_matrix_elements[i]);
    }

    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::LIE::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationLocalAssemblerMatrix,
        SmallDeformationLocalAssemblerMatrixNearFracture,
        SmallDeformationLocalAssemblerFracture>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    _secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXX));

    _secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYY));

    _secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaZZ));

    _secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        _secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXZ));

        _secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYZ));
    }

    _secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXX));

    _secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonYY));

    _secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonZZ));

    _secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXY));

    if (DisplacementDim == 3)
    {
        _secondary_variables.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXZ));

        _secondary_variables.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonYZ));
    }

    auto mesh_prop_sigma_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xx->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_stress_xx = mesh_prop_sigma_xx;

    auto mesh_prop_sigma_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_yy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_stress_yy = mesh_prop_sigma_yy;

    auto mesh_prop_sigma_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_zz->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_stress_zz = mesh_prop_sigma_zz;

    auto mesh_prop_sigma_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_stress_xy = mesh_prop_sigma_xy;

    if (DisplacementDim == 3)
    {
        auto mesh_prop_sigma_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_xz->resize(mesh.getNumberOfElements());
        _process_data._mesh_prop_stress_xz = mesh_prop_sigma_xz;

        auto mesh_prop_sigma_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_yz->resize(mesh.getNumberOfElements());
        _process_data._mesh_prop_stress_yz = mesh_prop_sigma_yz;
    }

    auto mesh_prop_epsilon_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xx->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_strain_xx = mesh_prop_epsilon_xx;

    auto mesh_prop_epsilon_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_yy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_strain_yy = mesh_prop_epsilon_yy;

    auto mesh_prop_epsilon_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_zz->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_strain_zz = mesh_prop_epsilon_zz;

    auto mesh_prop_epsilon_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_strain_xy = mesh_prop_epsilon_xy;

    if (DisplacementDim == 3)
    {
        auto mesh_prop_epsilon_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_xz->resize(mesh.getNumberOfElements());
        _process_data._mesh_prop_strain_xz = mesh_prop_epsilon_xz;

        auto mesh_prop_epsilon_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_yz->resize(mesh.getNumberOfElements());
        _process_data._mesh_prop_strain_yz = mesh_prop_epsilon_yz;
    }

    for (MeshLib::Element const* e : _mesh.getElements())
    {
        if (e->getDimension() < DisplacementDim)
            continue;

        Eigen::Vector3d const pt(e->getCenterOfGravity().getCoords());
        std::vector<FractureProperty*> e_fracture_props;
        std::unordered_map<int,int> e_fracID_to_local;
        unsigned tmpi = 0;
        for (auto fid : _process_data._vec_ele_connected_fractureIDs[e->getID()])
        {
            e_fracture_props.push_back(_process_data
                    ._vec_fracture_property[fid].get());
            e_fracID_to_local.insert({fid, tmpi++});
        }
        std::vector<JunctionProperty*> e_junction_props;
        std::unordered_map<int,int> e_juncID_to_local;
        tmpi = 0;
        for (auto fid : _process_data._vec_ele_connected_junctionIDs[e->getID()])
        {
            e_junction_props.push_back(
                _process_data._vec_junction_property[fid].get());
            e_juncID_to_local.insert({fid, tmpi++});
        }
        std::vector<double> const levelsets(
            u_global_enrichments(e_fracture_props, e_junction_props,
                                 e_fracID_to_local, pt));

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
                                   _process_data._vec_fracture_property.size()),
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_levelset->resize(mesh.getNumberOfElements());
            (*mesh_prop_levelset)[e->getID()] =
                levelsets[i + e_fracture_props.size()];
        }
    }


    auto mesh_prop_w_n = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "w_n", MeshLib::MeshItemType::Cell,
        1);
    mesh_prop_w_n->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_w_n = mesh_prop_w_n;

    auto mesh_prop_w_s = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "w_s", MeshLib::MeshItemType::Cell,
        1);
    mesh_prop_w_s->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_w_s = mesh_prop_w_s;

    auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "aperture",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_b->resize(mesh.getNumberOfElements());
    auto const& mesh_prop_matid = *_process_data._mesh_prop_materialIDs;
    for (auto const& fracture_prop : _process_data._vec_fracture_property)
    {
        for (MeshLib::Element const* e : _mesh.getElements())
        {
            if (e->getDimension() == DisplacementDim)
            {
                continue;
            }
            if (mesh_prop_matid[e->getID()] != fracture_prop->mat_id)
            {
                continue;
            }
            ProcessLib::SpatialPosition x;
            x.setElementID(e->getID());
            (*mesh_prop_b)[e->getID()] = (*fracture_prop->aperture0)(0, x)[0];
        }
    }
    _process_data._mesh_prop_b = mesh_prop_b;

    auto mesh_prop_fracture_stress_shear =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "f_stress_s",
            MeshLib::MeshItemType::Cell, 1);
    mesh_prop_fracture_stress_shear->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_fracture_stress_shear =
        mesh_prop_fracture_stress_shear;

    auto mesh_prop_fracture_stress_normal =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "f_stress_n",
            MeshLib::MeshItemType::Cell, 1);
    mesh_prop_fracture_stress_normal->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_fracture_stress_normal =
        mesh_prop_fracture_stress_normal;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x, int const process_id)
{
    DBUG("Compute the secondary variables for SmallDeformationProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &SmallDeformationLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, getDOFTable(process_id), t, x,
        _coupled_solutions);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble SmallDeformationProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}
template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    // Call global assembler for each local assembly item.
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_table, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, _coupled_solutions);
}
template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int /*process_id*/)
{
    DBUG("PreTimestep SmallDeformationProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &SmallDeformationLocalAssemblerInterface::preTimestep,
        _local_assemblers, *_local_to_global_index_map, x, t, dt);
}
// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class SmallDeformationProcess<2>;
template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
