/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationProcess-fwd.h"
#include "SmallDeformationProcess.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "ProcessLib/SmallDeformationWithLIE/BoundaryCondition/BoundaryConditionBuilder.h"
#include "ProcessLib/SmallDeformationWithLIE/Common/LevelSetFunction.h"
#include "ProcessLib/SmallDeformationWithLIE/Common/MeshUtils.h"
#include "ProcessLib/SmallDeformationWithLIE/Common/Utils.h"
#include "ProcessLib/SmallDeformationWithLIE/LocalAssembler/LocalAssemblerDataMatrix.h"
#include "ProcessLib/SmallDeformationWithLIE/LocalAssembler/LocalAssemblerDataMatrixNearFracture.h"
#include "ProcessLib/SmallDeformationWithLIE/LocalAssembler/LocalAssemblerDataFracture.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

template <int DisplacementDim>
SmallDeformationProcess<DisplacementDim>::SmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
        jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&&
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
    getFractureMatrixDataInMesh(mesh,
                                _vec_matrix_elements,
                                _vec_fracture_elements,
                                _vec_fracture_matrix_elements,
                                _vec_fracture_nodes);

    // set fracture property assuming a fracture forms a straight line
    setFractureProperty(DisplacementDim,
                        *_vec_fracture_elements[0],
                        *_process_data._fracture_property.get());

    // need to use a custom Neumann BC assembler for displacement jumps
    for (ProcessVariable& pv : getProcessVariables())
    {
        if (pv.getName().find("displacement_jump") == std::string::npos)
            continue;
        pv.setBoundaryConditionBuilder(
                    std::unique_ptr<ProcessLib::BoundaryConditionBuilder>(
                        new BoundaryConditionBuilder(*_process_data._fracture_property.get())));
    }

}


template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::constructDofTable()
{
    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // for extrapolation
    _mesh_subset_all_nodes.reset(new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));
    // regular u
    _mesh_subset_matrix_nodes.reset(new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));
    // u jump
    _mesh_subset_fracture_nodes.reset(new MeshLib::MeshSubset(_mesh, &_vec_fracture_nodes));

    // Collect the mesh subsets in a vector.
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
    std::generate_n(
        std::back_inserter(all_mesh_subsets),
        DisplacementDim,
        [&]() {
            return std::unique_ptr<MeshLib::MeshSubsets>{
                new MeshLib::MeshSubsets{_mesh_subset_matrix_nodes.get()}};
        });
    std::generate_n(
        std::back_inserter(all_mesh_subsets),
        DisplacementDim,
        [&]() {
            return std::unique_ptr<MeshLib::MeshSubsets>{
                new MeshLib::MeshSubsets{_mesh_subset_fracture_nodes.get()}};
        });

    std::vector<unsigned> const vec_n_components(2, DisplacementDim);

    std::vector<std::vector<MeshLib::Element*>const*> vec_var_elements;
    vec_var_elements.push_back(&_vec_matrix_elements);
    vec_var_elements.push_back(&_vec_fracture_matrix_elements);

    _local_to_global_index_map.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT));
}


template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformationWithLIE::createLocalAssemblers
            <DisplacementDim,
             LocalAssemblerDataMatrix,
             LocalAssemblerDataMatrixNearFracture,
             LocalAssemblerDataFracture>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>>
        all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        new MeshLib::MeshSubsets(_mesh_subset_all_nodes.get()));
    _local_to_global_index_map_single_component.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy", 1,
        makeExtrapolator(
            getExtrapolator(), _local_assemblers,
            &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXY));

    auto mesh_prop_sigma_xx = const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("stress_xx", MeshLib::MeshItemType::Cell);
    mesh_prop_sigma_xx->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_stress_xx = mesh_prop_sigma_xx;

    auto mesh_prop_sigma_yy = const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("stress_yy", MeshLib::MeshItemType::Cell);
    mesh_prop_sigma_yy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_stress_yy = mesh_prop_sigma_yy;

    auto mesh_prop_sigma_xy = const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("stress_xy", MeshLib::MeshItemType::Cell);
    mesh_prop_sigma_xy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_stress_xy = mesh_prop_sigma_xy;

    auto mesh_prop_epsilon_xx = const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("strain_xx", MeshLib::MeshItemType::Cell);
    mesh_prop_epsilon_xx->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_strain_xx = mesh_prop_epsilon_xx;

    auto mesh_prop_epsilon_yy = const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("strain_yy", MeshLib::MeshItemType::Cell);
    mesh_prop_epsilon_yy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_strain_yy = mesh_prop_epsilon_yy;

    auto mesh_prop_epsilon_xy = const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("strain_xy", MeshLib::MeshItemType::Cell);
    mesh_prop_epsilon_xy->resize(mesh.getNumberOfElements());
    _process_data._mesh_prop_strain_xy = mesh_prop_epsilon_xy;

    auto mesh_prop_levelset =
            const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("levelset1", MeshLib::MeshItemType::Cell);
    mesh_prop_levelset->resize(mesh.getNumberOfElements());
    for (MeshLib::Element const* e : _mesh.getElements())
    {
        if (e->getDimension() < DisplacementDim)
            continue;

        double const levelsets = calculateLevelSetFunction(*_process_data._fracture_property, e->getCenterOfGravity().getCoords());
        (*mesh_prop_levelset)[e->getID()] = levelsets;
    }

    auto mesh_prop_b = const_cast<MeshLib::Mesh&>(mesh).getProperties().template createNewPropertyVector<double>("aperture", MeshLib::MeshItemType::Cell);
    mesh_prop_b->resize(mesh.getNumberOfElements());
    auto mesh_prop_matid = mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    auto frac = _process_data._fracture_property.get();
    for (MeshLib::Element const* e : _mesh.getElements())
    {
        if (e->getDimension() == DisplacementDim)
            continue;
        if ((*mesh_prop_matid)[e->getID()] != frac->mat_id)
            continue;
        ProcessLib::SpatialPosition x;
        x.setElementID(e->getID());
        (*mesh_prop_b)[e->getID()] = (*frac->aperture0)(0, x)[0];
    }
    _process_data._mesh_prop_b = mesh_prop_b;
}


template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::postTimestepConcreteProcess(GlobalVector const& x)
{
    DBUG("PostTimestep SmallDeformationProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &SmallDeformationLocalAssemblerInterface::postTimestep,
        _local_assemblers, *_local_to_global_index_map, x);

}


// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class SmallDeformationProcess<2>;

}   // namespace SmallDeformationWithLIE
}   // namespace ProcessLib
