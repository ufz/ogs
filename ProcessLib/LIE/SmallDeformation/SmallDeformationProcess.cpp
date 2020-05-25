/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      process_data_(std::move(process_data))
{
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_branch_nodeID_matIDs;
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_junction_nodeID_matIDs;
    getFractureMatrixDataInMesh(mesh, vec_matrix_elements_,
                                vec_fracture_mat_IDs_, vec_fracture_elements_,
                                vec_fracture_matrix_elements_,
                                vec_fracture_nodes_, vec_branch_nodeID_matIDs,
                                vec_junction_nodeID_matIDs);

    if (vec_fracture_mat_IDs_.size() !=
        process_data_.fracture_properties.size())
    {
        OGS_FATAL(
            "The number of the given fracture properties ({:d}) are not "
            "consistent"
            " with the number of fracture groups in a mesh ({:d}).",
            process_data_.fracture_properties.size(),
            vec_fracture_mat_IDs_.size());
    }

    // create a map from a material ID to a fracture ID
    auto max_frac_mat_id = std::max_element(vec_fracture_mat_IDs_.begin(),
                                            vec_fracture_mat_IDs_.end());
    process_data_.map_materialID_to_fractureID_.resize(*max_frac_mat_id + 1);
    for (unsigned i = 0; i < vec_fracture_mat_IDs_.size(); i++)
    {
        process_data_.map_materialID_to_fractureID_[vec_fracture_mat_IDs_[i]] =
            i;
    }

    // create a table of connected fracture IDs for each element
    process_data_.vec_ele_connected_fractureIDs_.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < vec_fracture_matrix_elements_.size(); i++)
    {
        for (auto e : vec_fracture_matrix_elements_[i])
        {
            process_data_.vec_ele_connected_fractureIDs_[e->getID()].push_back(
                i);
        }
    }

    // set fracture property
    for (auto& fracture_prop : process_data_.fracture_properties)
    {
        // based on the 1st element assuming a fracture forms a straight line
        setFractureProperty(
            DisplacementDim,
            *vec_fracture_elements_[fracture_prop.fracture_id][0],
            fracture_prop);
    }

    // set branches
    for (auto& vec_branch_nodeID_matID : vec_branch_nodeID_matIDs)
    {
        auto master_matId = vec_branch_nodeID_matID.second[0];
        auto slave_matId = vec_branch_nodeID_matID.second[1];
        auto& master_frac =
            process_data_.fracture_properties
                [process_data_.map_materialID_to_fractureID_[master_matId]];
        auto& slave_frac =
            process_data_.fracture_properties
                [process_data_.map_materialID_to_fractureID_[slave_matId]];

        master_frac.branches_master.push_back(
            createBranchProperty(*mesh.getNode(vec_branch_nodeID_matID.first),
                                 master_frac, slave_frac));

        slave_frac.branches_slave.push_back(
            createBranchProperty(*mesh.getNode(vec_branch_nodeID_matID.first),
                                 master_frac, slave_frac));
    }

    // set junctions
    for (auto& vec_junction_nodeID_matID : vec_junction_nodeID_matIDs)
    {
        vec_junction_nodes_.push_back(const_cast<MeshLib::Node*>(
            mesh_.getNode(vec_junction_nodeID_matID.first)));
    }
    for (std::size_t i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto const& material_ids = vec_junction_nodeID_matIDs[i].second;
        assert(material_ids.size() == 2);
        std::array<int, 2> fracture_ids{
            {process_data_.map_materialID_to_fractureID_[material_ids[0]],
             process_data_.map_materialID_to_fractureID_[material_ids[1]]}};

        process_data_.junction_properties.emplace_back(
            i, *mesh.getNode(vec_junction_nodeID_matIDs[i].first),
            fracture_ids);
    }

    // create a table of connected junction IDs for each element
    process_data_.vec_ele_connected_junctionIDs_.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto e : node->getElements())
        {
            process_data_.vec_ele_connected_junctionIDs_[e->getID()].push_back(
                i);
        }
    }

    // create a table of junction node and connected elements
    vec_junction_fracture_matrix_elements_.resize(
        vec_junction_nodeID_matIDs.size());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto e : node->getElements())
        {
            vec_junction_fracture_matrix_elements_[i].push_back(e);
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
    process_data_.mesh_prop_materialIDs_ = material_ids;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::constructDofTable()
{
    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // for extrapolation
    mesh_subset_all_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, mesh_.getNodes());
    // regular u
    mesh_subset_matrix_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, mesh_.getNodes());
    // u jump
    for (unsigned i = 0; i < vec_fracture_nodes_.size(); i++)
    {
        mesh_subset_fracture_nodes_.push_back(
            std::make_unique<MeshLib::MeshSubset const>(
                mesh_, vec_fracture_nodes_[i]));
    }
    // enrichment for junctions
    mesh_subset_junction_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, vec_junction_nodes_);

    // Collect the mesh subsets in a vector.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::generate_n(std::back_inserter(all_mesh_subsets), DisplacementDim,
                    [&]() { return *mesh_subset_matrix_nodes_; });
    for (auto& ms : mesh_subset_fracture_nodes_)
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        DisplacementDim,
                        [&]() { return *ms; });
    }
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    DisplacementDim,
                    [&]() { return *mesh_subset_junction_nodes_; });

    std::vector<int> const vec_n_components(
        1 + vec_fracture_mat_IDs_.size() + vec_junction_nodes_.size(),
        DisplacementDim);

    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    vec_var_elements.push_back(&vec_matrix_elements_);
    for (unsigned i = 0; i < vec_fracture_matrix_elements_.size(); i++)
    {
        vec_var_elements.push_back(&vec_fracture_matrix_elements_[i]);
    }
    for (unsigned i = 0; i < vec_junction_fracture_matrix_elements_.size(); i++)
    {
        vec_var_elements.push_back(&vec_junction_fracture_matrix_elements_[i]);
    }

    local_to_global_index_map_ =
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
        mesh.getElements(), dof_table, local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, process_data_);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *mesh_subset_all_nodes_};
    local_to_global_index_map_single_component_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    secondary_variables_.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXX));

    secondary_variables_.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYY));

    secondary_variables_.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaZZ));

    secondary_variables_.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        secondary_variables_.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(
                1, getExtrapolator(), local_assemblers_,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXZ));

        secondary_variables_.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(
                1, getExtrapolator(), local_assemblers_,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYZ));
    }

    secondary_variables_.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXX));

    secondary_variables_.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonYY));

    secondary_variables_.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonZZ));

    secondary_variables_.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(
            1, getExtrapolator(), local_assemblers_,
            &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXY));

    if (DisplacementDim == 3)
    {
        secondary_variables_.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(
                1, getExtrapolator(), local_assemblers_,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXZ));

        secondary_variables_.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(
                1, getExtrapolator(), local_assemblers_,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonYZ));
    }

    auto mesh_prop_sigma_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xx->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_stress_xx_ = mesh_prop_sigma_xx;

    auto mesh_prop_sigma_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_yy->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_stress_yy_ = mesh_prop_sigma_yy;

    auto mesh_prop_sigma_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_zz->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_stress_zz_ = mesh_prop_sigma_zz;

    auto mesh_prop_sigma_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xy->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_stress_xy_ = mesh_prop_sigma_xy;

    if (DisplacementDim == 3)
    {
        auto mesh_prop_sigma_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_xz->resize(mesh.getNumberOfElements());
        process_data_.mesh_prop_stress_xz_ = mesh_prop_sigma_xz;

        auto mesh_prop_sigma_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_yz->resize(mesh.getNumberOfElements());
        process_data_.mesh_prop_stress_yz_ = mesh_prop_sigma_yz;
    }

    auto mesh_prop_epsilon_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xx->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_strain_xx_ = mesh_prop_epsilon_xx;

    auto mesh_prop_epsilon_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_yy->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_strain_yy_ = mesh_prop_epsilon_yy;

    auto mesh_prop_epsilon_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_zz->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_strain_zz_ = mesh_prop_epsilon_zz;

    auto mesh_prop_epsilon_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xy->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_strain_xy_ = mesh_prop_epsilon_xy;

    if (DisplacementDim == 3)
    {
        auto mesh_prop_epsilon_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_xz->resize(mesh.getNumberOfElements());
        process_data_.mesh_prop_strain_xz_ = mesh_prop_epsilon_xz;

        auto mesh_prop_epsilon_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_yz->resize(mesh.getNumberOfElements());
        process_data_.mesh_prop_strain_yz_ = mesh_prop_epsilon_yz;
    }

    for (MeshLib::Element const* e : mesh_.getElements())
    {
        if (e->getDimension() < DisplacementDim)
        {
            continue;
        }

        Eigen::Vector3d const pt(e->getCenterOfGravity().getCoords());
        std::vector<FractureProperty*> e_fracture_props;
        std::unordered_map<int, int> e_fracID_to_local;
        unsigned tmpi = 0;
        for (auto fid :
             process_data_.vec_ele_connected_fractureIDs_[e->getID()])
        {
            e_fracture_props.push_back(&process_data_.fracture_properties[fid]);
            e_fracID_to_local.insert({fid, tmpi++});
        }
        std::vector<JunctionProperty*> e_junction_props;
        std::unordered_map<int, int> e_juncID_to_local;
        tmpi = 0;
        for (auto fid :
             process_data_.vec_ele_connected_junctionIDs_[e->getID()])
        {
            e_junction_props.push_back(&process_data_.junction_properties[fid]);
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
                                   process_data_.fracture_properties.size()),
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
    process_data_.mesh_prop_w_n_ = mesh_prop_w_n;

    auto mesh_prop_w_s = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "w_s", MeshLib::MeshItemType::Cell,
        1);
    mesh_prop_w_s->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_w_s_ = mesh_prop_w_s;

    auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "aperture",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_b->resize(mesh.getNumberOfElements());
    auto const& mesh_prop_matid = *process_data_.mesh_prop_materialIDs_;
    for (auto const& fracture_prop : process_data_.fracture_properties)
    {
        for (MeshLib::Element const* e : mesh_.getElements())
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
    process_data_.mesh_prop_b_ = mesh_prop_b;

    auto mesh_prop_fracture_stress_shear =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "f_stress_s",
            MeshLib::MeshItemType::Cell, 1);
    mesh_prop_fracture_stress_shear->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_fracture_stress_shear_ =
        mesh_prop_fracture_stress_shear;

    auto mesh_prop_fracture_stress_normal =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "f_stress_n",
            MeshLib::MeshItemType::Cell, 1);
    mesh_prop_fracture_stress_normal->resize(mesh.getNumberOfElements());
    process_data_.mesh_prop_fracture_stress_normal_ =
        mesh_prop_fracture_stress_normal;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, GlobalVector const& x,
    GlobalVector const& x_dot, int const process_id)
{
    DBUG("Compute the secondary variables for SmallDeformationProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &SmallDeformationLocalAssemblerInterface::computeSecondaryVariable,
        local_assemblers_, pv.getActiveElementIDs(), getDOFTable(process_id), t,
        dt, x, x_dot, coupled_solutions_);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble SmallDeformationProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}
template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    // Call global assembler for each local assembly item.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*local_to_global_index_map_)};
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);
}
template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep SmallDeformationProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &SmallDeformationLocalAssemblerInterface::preTimestep,
        local_assemblers_, pv.getActiveElementIDs(),
        *local_to_global_index_map_, *x[process_id], t, dt);
}
// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class SmallDeformationProcess<2>;
template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
