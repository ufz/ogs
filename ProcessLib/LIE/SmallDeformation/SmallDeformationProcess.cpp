/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationProcess.h"

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
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
      _process_data(std::move(process_data))
{
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
                [_process_data._map_materialID_to_fractureID[master_matId]];
        auto& slave_frac =
            _process_data.fracture_properties
                [_process_data._map_materialID_to_fractureID[slave_matId]];

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
            {_process_data._map_materialID_to_fractureID[material_ids[0]],
             _process_data._map_materialID_to_fractureID[material_ids[1]]}};

        _process_data.junction_properties.emplace_back(
            i, *mesh.getNode(vec_junction_nodeID_matIDs[i].first),
            fracture_ids);
    }

    // create a table of connected junction IDs for each element
    _process_data._vec_ele_connected_junctionIDs.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto id :
             mesh.getElementsConnectedToNode(*node) | MeshLib::views::ids)
        {
            _process_data._vec_ele_connected_junctionIDs[id].push_back(i);
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
    for (auto const& ms : _mesh_subset_fracture_nodes)
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        DisplacementDim,
                        [&]() { return *ms; });
    }
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    DisplacementDim,
                    [&]() { return *_mesh_subset_junction_nodes; });

    std::vector<int> const vec_n_components(
        1 + _vec_fracture_mat_IDs.size() + _vec_junction_nodes.size(),
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
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

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

    add_secondary_variable("fracture_stress", DisplacementDim,
                           &LocalAssemblerInterface::getIntPtFractureStress);

    add_secondary_variable("fracture_aperture", 1,
                           &LocalAssemblerInterface::getIntPtFractureAperture);

    _process_data.element_stresses = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "sigma_avg",
        MeshLib::MeshItemType::Cell,
        MathLib::KelvinVector::KelvinVectorType<
            DisplacementDim>::RowsAtCompileTime);

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
        for (auto fid :
             _process_data._vec_ele_connected_fractureIDs[e->getID()])
        {
            e_fracture_props.push_back(&_process_data.fracture_properties[fid]);
            e_fracID_to_local.insert({fid, tmpi++});
        }
        std::vector<JunctionProperty*> e_junction_props;
        std::unordered_map<int, int> e_juncID_to_local;
        tmpi = 0;
        for (auto fid :
             _process_data._vec_ele_connected_junctionIDs[e->getID()])
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

    auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "fracture_aperture_avg",
        MeshLib::MeshItemType::Cell, 1);

    mesh_prop_b->resize(mesh.getNumberOfElements());
    auto const& mesh_prop_matid = *_process_data._mesh_prop_materialIDs;
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
    _process_data._mesh_prop_b = mesh_prop_b;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, int const process_id)
{
    DBUG("Compute the secondary variables for SmallDeformationProcess.");
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &SmallDeformationLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_prev, process_id);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble SmallDeformationProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, M, K,
        b);
}
template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    // Call global assembler for each local assembly item.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x,
        x_prev, process_id, M, K, b, Jac);
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
        _local_assemblers, pv.getActiveElementIDs(),
        *_local_to_global_index_map, *x[process_id], t, dt);
}
// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class SmallDeformationProcess<2>;
template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
