/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonBoundaryCondition.h"

#include <pybind11/pybind11.h>

#include <iostream>

#include "BaseLib/ConfigTree.h"
#include "FlushStdoutGuard.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/ProcessVariable.h"
#include "PythonBoundaryConditionLocalAssembler.h"

namespace
{
void initBCValues(NumLib::IndexValueVector<GlobalIndexType>& bc_values,
                  std::size_t const nnodes)
{
    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.reserve(nnodes);
    bc_values.values.reserve(nnodes);
}

void checkConsistency(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>> const& pvs)
{
    auto const num_vars_dt = dof_table.getNumberOfVariables();
    auto const num_vars_pv = pvs.size();

    if (static_cast<std::size_t>(num_vars_dt) != num_vars_pv)
    {
        OGS_FATAL(
            "The number of variables in the d.o.f. table does not match the "
            "number of process variables: {} != {}.",
            num_vars_dt, num_vars_pv);
    }

    for (std::size_t var = 0; var < num_vars_pv; ++var)
    {
        auto const num_comp_dt = dof_table.getNumberOfVariableComponents(var);
        auto const num_comp_pv = pvs[var].get().getNumberOfGlobalComponents();

        if (num_comp_dt != num_comp_pv)
        {
            OGS_FATAL(
                "The number of components of variable #{} in the d.o.f. table "
                "does not match the number of components of process variable "
                "#{} ({}): {} != {}.",
                var, var, pvs[var].get().getName(), num_vars_dt, num_vars_pv);
        }
    }
}
}  // anonymous namespace

namespace ProcessLib
{
PythonBoundaryCondition::PythonBoundaryCondition(
    PythonBcData&& bc_data, unsigned const integration_order,
    bool const flush_stdout, unsigned const bulk_mesh_dimension,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk)
    : _bc_data(std::move(bc_data)), _flush_stdout(flush_stdout)
{
    checkConsistency(dof_table_bulk,
                     _bc_data.all_process_variables_for_this_process);

    std::vector<MeshLib::Node*> const& bc_nodes =
        _bc_data.bc_or_st_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(_bc_data.bc_or_st_mesh, bc_nodes);

    _dof_table_boundary =
        dof_table_bulk.deriveBoundaryConstrainedMap(std::move(bc_mesh_subset));

    BoundaryConditionAndSourceTerm::createLocalAssemblersPython<
        PythonBoundaryConditionLocalAssembler>(
        bulk_mesh_dimension, _bc_data.bc_or_st_mesh.getElements(),
        *_dof_table_boundary, _local_assemblers,
        _bc_data.bc_or_st_mesh.isAxiallySymmetric(), integration_order,
        _bc_data);
}

void PythonBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    FlushStdoutGuard guard(_flush_stdout);
    (void)guard;

    auto const& boundary_mesh = _bc_data.bc_or_st_mesh;
    auto const& boundary_nodes = boundary_mesh.getNodes();
    auto const* bc_object = _bc_data.bc_or_st_object;

    initBCValues(bc_values, boundary_nodes.size());

    std::vector<double> primary_variables;

    for (auto const* boundary_node : boundary_nodes)
    {
        auto const boundary_node_id = boundary_node->getID();
        auto const dof_idx = getDofIdx(boundary_node_id);

        if (dof_idx == NumLib::MeshComponentMap::nop)
        {
            // This d.o.f. has not been found. This can be the case, e.g., for
            // Taylor-Hood elements, where the lower order field has no d.o.f.
            // on higher order nodes.
            continue;
        }

        if (dof_idx < 0)
        {
            // For the DDC approach (e.g. with PETSc option) a negative
            // index means that this entry is a ghost entry and should be
            // dropped.
            continue;
        }

        collectPrimaryVariables(primary_variables, *boundary_node, x);

        auto const [apply_bc, bc_value] = bc_object->getDirichletBCValue(
            t, {(*boundary_node)[0], (*boundary_node)[1], (*boundary_node)[2]},
            boundary_node_id, primary_variables);

        if (!bc_object->isOverriddenEssential())
        {
            DBUG(
                "Method `getDirichletBCValue' not overridden in Python "
                "script.");
            return;
        }

        if (!apply_bc)
        {
            continue;
        }

        bc_values.ids.emplace_back(dof_idx);
        bc_values.values.emplace_back(bc_value);
    }
}

GlobalIndexType PythonBoundaryCondition::getDofIdx(
    std::size_t const boundary_node_id) const
{
    MeshLib::Location const loc{_bc_data.bc_or_st_mesh.getID(),
                                MeshLib::MeshItemType::Node, boundary_node_id};
    return _dof_table_boundary->getGlobalIndex(loc,
                                               _bc_data.global_component_id);
}

GlobalIndexType PythonBoundaryCondition::getDofIdx(
    std::size_t const boundary_node_id, int const var, int const comp) const
{
    MeshLib::Location const loc{_bc_data.bc_or_st_mesh.getID(),
                                MeshLib::MeshItemType::Node, boundary_node_id};
    return _dof_table_boundary->getGlobalIndex(loc, var, comp);
}

void PythonBoundaryCondition::collectPrimaryVariables(
    std::vector<double>& primary_variables, MeshLib::Node const& boundary_node,
    GlobalVector const& x) const
{
    primary_variables.clear();
    auto const num_var = _dof_table_boundary->getNumberOfVariables();
    auto const boundary_node_id = boundary_node.getID();

    for (int var = 0; var < num_var; ++var)
    {
        auto const num_comp =
            _dof_table_boundary->getNumberOfVariableComponents(var);
        for (int comp = 0; comp < num_comp; ++comp)
        {
            auto const dof_idx = getDofIdx(boundary_node_id, var, comp);

            double const pv_value =
                dof_idx != NumLib::MeshComponentMap::nop
                    ? x[dof_idx]
                    : interpolateToHigherOrderNode(x, var, comp, boundary_node);

            primary_variables.push_back(pv_value);
        }
    }
}

double PythonBoundaryCondition::interpolateToHigherOrderNode(
    GlobalVector const& x, int const var, int const comp,
    MeshLib::Node const& boundary_node) const
{
    auto const& boundary_elements =
        _bc_data.bc_or_st_mesh.getElementsConnectedToNode(
            boundary_node.getID());

    if (boundary_elements.size() != 1)
    {
        DBUG(
            "Boundary node {} is associated with {} elements in the boundary "
            "mesh.",
            boundary_node.getID(),
            boundary_elements.size());
    }

    // Interpolations on all associated elements should return the same. Just
    // pick any of them.
    auto const& boundary_element = *boundary_elements.front();

    assert(boundary_element.getNumberOfBaseNodes() <
               boundary_element.getNumberOfNodes() &&
           "We expect that the boundary element is a higher order element. "
           "Otherwise no interpolation should take place.");

    // Search local node id
    auto const local_node_id_within_boundary_element =
        getNodeIDinElement(boundary_element, &boundary_node);

    auto const boundary_element_id = boundary_element.getID();

    // Assumption: all boundary elements have a local assembler (in the same
    // order)
    auto const& loc_asm = *_local_assemblers[boundary_element_id];

    return loc_asm.interpolate(local_node_id_within_boundary_element,
                               *_dof_table_boundary, x, var, comp);
}

void PythonBoundaryCondition::applyNaturalBC(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    GlobalMatrix& K, GlobalVector& b, GlobalMatrix* Jac)
{
    FlushStdoutGuard guard(_flush_stdout);

    try
    {
        GlobalExecutor::executeMemberOnDereferenced(
            &GenericNaturalBoundaryConditionLocalAssemblerInterface::assemble,
            _local_assemblers, *_dof_table_boundary, t, x, process_id, K, b,
            Jac);
    }
    catch (MethodNotOverriddenInDerivedClassException const& /*e*/)
    {
        DBUG("Method `getFlux' not overridden in Python script.");
    }
}

std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    MeshLib::Mesh const& bulk_mesh, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process)
{
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "Python");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Python__bc_object}
    auto const bc_object = config.getConfigParameter<std::string>("bc_object");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Python__flush_stdout}
    auto const flush_stdout = config.getConfigParameter("flush_stdout", false);

    // Evaluate Python code in scope of main module
    pybind11::object scope =
        pybind11::module::import("__main__").attr("__dict__");

    if (!scope.contains(bc_object))
    {
        OGS_FATAL(
            "Function `{:s}' is not defined in the python script file, or "
            "there was no python script file specified.",
            bc_object);
    }

    auto* bc = scope[bc_object.c_str()]
                   .cast<PythonBoundaryConditionPythonSideInterface*>();

    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: ({:d}, "
            "{:d}), maximum values: ({:d}, {:d}).",
            variable_id, component_id, dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    // In case of partitioned mesh the boundary could be empty, i.e. there is no
    // boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createBoundaryCondition().
    if (boundary_mesh.getDimension() == 0 &&
        boundary_mesh.getNumberOfNodes() == 0 &&
        boundary_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<PythonBoundaryCondition>(
        PythonBcData{
            {bc, dof_table_bulk.getGlobalComponent(variable_id, component_id),
             boundary_mesh, all_process_variables_for_this_process,
             shapefunction_order}},
        integration_order, flush_stdout, bulk_mesh.getDimension(),
        dof_table_bulk);
}

}  // namespace ProcessLib
