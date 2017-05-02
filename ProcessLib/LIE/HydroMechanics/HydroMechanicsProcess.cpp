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
            std::unique_ptr<ProcessLib::BoundaryConditionBuilder>(
                new BoundaryConditionBuilder(
                    *_process_data.fracture_property.get())));
    }

    if (!_process_data.deactivate_matrix_in_flow)
    {
        _process_data.p_element_status.reset(new MeshLib::ElementStatus(&mesh));
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
        _process_data.p_element_status.reset(new MeshLib::ElementStatus(&mesh, vec_p_inactive_matIDs));

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
    _mesh_subset_all_nodes.reset(
        new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));
    // pressure
    _mesh_nodes_p = MeshLib::getBaseNodes(
        _process_data.p_element_status->getActiveElements());
    _mesh_subset_nodes_p.reset(new MeshLib::MeshSubset(_mesh, &_mesh_nodes_p));
    // regular u
    _mesh_subset_matrix_nodes.reset(
        new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));
    if (!_vec_fracture_nodes.empty())
    {
        // u jump
        _mesh_subset_fracture_nodes.reset(
            new MeshLib::MeshSubset(_mesh, &_vec_fracture_nodes));
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
    _local_to_global_index_map.reset(new NumLib::LocalToGlobalIndexMap(
        std::move(all_mesh_subsets),
        vec_n_components,
        vec_var_elements,
        NumLib::ComponentOrder::BY_COMPONENT));

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

    auto mesh_prop_sigma_xx = const_cast<MeshLib::Mesh&>(mesh)
                                  .getProperties()
                                  .template createNewPropertyVector<double>(
                                      "stress_xx", MeshLib::MeshItemType::Cell);
    mesh_prop_sigma_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xx = mesh_prop_sigma_xx;

    auto mesh_prop_sigma_yy = const_cast<MeshLib::Mesh&>(mesh)
                                  .getProperties()
                                  .template createNewPropertyVector<double>(
                                      "stress_yy", MeshLib::MeshItemType::Cell);
    mesh_prop_sigma_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_yy = mesh_prop_sigma_yy;

    auto mesh_prop_sigma_xy = const_cast<MeshLib::Mesh&>(mesh)
                                  .getProperties()
                                  .template createNewPropertyVector<double>(
                                      "stress_xy", MeshLib::MeshItemType::Cell);
    mesh_prop_sigma_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xy = mesh_prop_sigma_xy;

    auto mesh_prop_epsilon_xx =
        const_cast<MeshLib::Mesh&>(mesh)
            .getProperties()
            .template createNewPropertyVector<double>(
                "strain_xx", MeshLib::MeshItemType::Cell);
    mesh_prop_epsilon_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xx = mesh_prop_epsilon_xx;

    auto mesh_prop_epsilon_yy =
        const_cast<MeshLib::Mesh&>(mesh)
            .getProperties()
            .template createNewPropertyVector<double>(
                "strain_yy", MeshLib::MeshItemType::Cell);
    mesh_prop_epsilon_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_yy = mesh_prop_epsilon_yy;

    auto mesh_prop_epsilon_xy =
        const_cast<MeshLib::Mesh&>(mesh)
            .getProperties()
            .template createNewPropertyVector<double>(
                "strain_xy", MeshLib::MeshItemType::Cell);
    mesh_prop_epsilon_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xy = mesh_prop_epsilon_xy;

    auto mesh_prop_velocity =
        const_cast<MeshLib::Mesh&>(mesh)
            .getProperties()
            .template createNewPropertyVector<double>(
                "velocity", MeshLib::MeshItemType::Cell, 3);
    mesh_prop_velocity->resize(mesh.getNumberOfElements() * 3);
    _process_data.mesh_prop_velocity = mesh_prop_velocity;

    if (!_vec_fracture_elements.empty())
    {
        auto mesh_prop_levelset =
            const_cast<MeshLib::Mesh&>(mesh)
                .getProperties()
                .template createNewPropertyVector<double>(
                    "levelset1", MeshLib::MeshItemType::Cell);
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

        auto mesh_prop_w_n = const_cast<MeshLib::Mesh&>(mesh)
                                 .getProperties()
                                 .template createNewPropertyVector<double>(
                                     "w_n", MeshLib::MeshItemType::Cell);
        mesh_prop_w_n->resize(mesh.getNumberOfElements());
        auto mesh_prop_w_s = const_cast<MeshLib::Mesh&>(mesh)
                                 .getProperties()
                                 .template createNewPropertyVector<double>(
                                     "w_s", MeshLib::MeshItemType::Cell);
        mesh_prop_w_s->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_w_n = mesh_prop_w_n;
        _process_data.mesh_prop_w_s = mesh_prop_w_s;

        auto mesh_prop_b = const_cast<MeshLib::Mesh&>(mesh)
                               .getProperties()
                               .template createNewPropertyVector<double>(
                                   "aperture", MeshLib::MeshItemType::Cell);
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

        auto mesh_prop_k_f = const_cast<MeshLib::Mesh&>(mesh)
                                 .getProperties()
                                 .template createNewPropertyVector<double>(
                                     "k_f", MeshLib::MeshItemType::Cell);
        mesh_prop_k_f->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_k_f = mesh_prop_k_f;

        auto mesh_prop_fracture_stress_shear =
            const_cast<MeshLib::Mesh&>(mesh)
                .getProperties()
                .template createNewPropertyVector<double>(
                    "f_stress_s", MeshLib::MeshItemType::Cell);
        mesh_prop_fracture_stress_shear->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_shear =
            mesh_prop_fracture_stress_shear;

        auto mesh_prop_fracture_stress_normal =
            const_cast<MeshLib::Mesh&>(mesh)
                .getProperties()
                .template createNewPropertyVector<double>(
                    "f_stress_n", MeshLib::MeshItemType::Cell);
        mesh_prop_fracture_stress_normal->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_normal =
            mesh_prop_fracture_stress_normal;

        auto mesh_prop_fracture_shear_failure =
            const_cast<MeshLib::Mesh&>(mesh)
                .getProperties()
                .template createNewPropertyVector<double>(
                    "f_shear_failure", MeshLib::MeshItemType::Cell);
        mesh_prop_fracture_shear_failure->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_shear_failure =
            mesh_prop_fracture_shear_failure;
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
}

// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class HydroMechanicsProcess<2>;

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
