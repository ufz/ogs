/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonBoundaryCondition.h"

#include <pybind11/pybind11.h>

#include <iostream>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "PythonBoundaryConditionLocalAssembler.h"

namespace
{
//! Optionally flushes the standard output upon creation and destruction.
//! Can be used to improve the debug output readability when printing debug
//! messages both from OGS and from Python.
class FlushStdoutGuard
{
public:
    //! Optionally flushes C++ stdout before running Python code.
    explicit FlushStdoutGuard(bool const flush) : flush_(flush)
    {
        if (!flush)
        {
            return;
        }

        std::cout << std::flush;
    }

    //! Optionally flushes Python's stdout after running Python code.
    ~FlushStdoutGuard()
    {
        if (!flush_)
        {
            return;
        }

        using namespace pybind11::literals;
        pybind11::print("end"_a = "", "flush"_a = true);
    }

private:
    //! To flush or not to flush.
    const bool flush_;
};
}  // anonymous namespace

namespace ProcessLib
{
PythonBoundaryCondition::PythonBoundaryCondition(
    PythonBoundaryConditionData&& bc_data,
    unsigned const integration_order,
    unsigned const shapefunction_order,
    unsigned const global_dim,
    bool const flush_stdout)
    : bc_data_(std::move(bc_data)), flush_stdout_(flush_stdout)
{
    std::vector<MeshLib::Node*> const& bc_nodes =
        bc_data_.boundary_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(bc_data_.boundary_mesh, bc_nodes);

    // Create local DOF table from the bc mesh subset for the given variable and
    // component id.
    dof_table_boundary_ = bc_data_.dof_table_bulk.deriveBoundaryConstrainedMap(
        std::move(bc_mesh_subset));

    createLocalAssemblers<PythonBoundaryConditionLocalAssembler>(
        global_dim, bc_data_.boundary_mesh.getElements(), *dof_table_boundary_,
        shapefunction_order, local_assemblers_,
        bc_data_.boundary_mesh.isAxiallySymmetric(), integration_order,
        bc_data_);
}

void PythonBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    FlushStdoutGuard guard(flush_stdout_);
    (void)guard;

    auto const nodes = bc_data_.boundary_mesh.getNodes();

    auto const& bulk_node_ids_map =
        *bc_data_.boundary_mesh.getProperties().getPropertyVector<std::size_t>(
            "bulk_node_ids", MeshLib::MeshItemType::Node, 1);

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.reserve(bc_data_.boundary_mesh.getNumberOfNodes());
    bc_values.values.reserve(bc_data_.boundary_mesh.getNumberOfNodes());

    std::vector<double> primary_variables;

    for (auto const* node : bc_data_.boundary_mesh.getNodes())
    {
        auto const boundary_node_id = node->getID();
        auto const bulk_node_id = bulk_node_ids_map[boundary_node_id];

        // gather primary variable values
        primary_variables.clear();
        auto const num_var = dof_table_boundary_->getNumberOfVariables();
        for (int var = 0; var < num_var; ++var)
        {
            auto const num_comp =
                dof_table_boundary_->getNumberOfVariableComponents(var);
            for (int comp = 0; comp < num_comp; ++comp)
            {
                MeshLib::Location loc{bc_data_.bulk_mesh_id,
                                      MeshLib::MeshItemType::Node,
                                      bulk_node_id};
                auto const dof_idx =
                    bc_data_.dof_table_bulk.getGlobalIndex(loc, var, comp);

                if (dof_idx == NumLib::MeshComponentMap::nop)
                {
                    // TODO extend Python BC to mixed FEM ansatz functions
                    OGS_FATAL(
                        "No d.o.f. found for (node={:d}, var={:d}, comp={:d}). "
                        " "
                        "That might be due to the use of mixed FEM ansatz "
                        "functions, which is currently not supported by "
                        "the implementation of Python BCs. That excludes, "
                        "e.g., the HM process.",
                        bulk_node_id, var, comp);
                }

                primary_variables.push_back(x[dof_idx]);
            }
        }

        auto* xs = nodes[boundary_node_id]->getCoords();  // TODO DDC problems?
        auto pair_flag_value = bc_data_.bc_object->getDirichletBCValue(
            t, {xs[0], xs[1], xs[2]}, boundary_node_id, primary_variables);
        if (!bc_data_.bc_object->isOverriddenEssential())
        {
            DBUG(
                "Method `getDirichletBCValue' not overridden in Python "
                "script.");
            return;
        }

        if (!pair_flag_value.first)
        {
            continue;
        }

        MeshLib::Location l(bc_data_.bulk_mesh_id, MeshLib::MeshItemType::Node,
                            bulk_node_id);
        const auto dof_idx = bc_data_.dof_table_bulk.getGlobalIndex(
            l, bc_data_.global_component_id);
        if (dof_idx == NumLib::MeshComponentMap::nop)
        {
            OGS_FATAL(
                "Logic error. This error should already have occured while "
                "gathering primary variables. Something nasty is going on!");
        }

        // For the DDC approach (e.g. with PETSc option), the negative
        // index of g_idx means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply
        // the Dirichlet BC, the negative index is not accepted like
        // other matrix or vector PETSc routines. Therefore, the
        // following if-condition is applied.
        if (dof_idx >= 0)
        {
            bc_values.ids.emplace_back(dof_idx);
            bc_values.values.emplace_back(pair_flag_value.second);
        }
    }
}

void PythonBoundaryCondition::applyNaturalBC(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    GlobalMatrix& K, GlobalVector& b, GlobalMatrix* Jac)
{
    FlushStdoutGuard guard(flush_stdout_);

    try
    {
        GlobalExecutor::executeMemberOnDereferenced(
            &GenericNaturalBoundaryConditionLocalAssemblerInterface::assemble,
            local_assemblers_, *dof_table_boundary_, t, x, process_id, K, b,
            Jac);
    }
    catch (MethodNotOverriddenInDerivedClassException const& /*e*/)
    {
        DBUG("Method `getFlux' not overridden in Python script.");
    }
}

std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t bulk_mesh_id,
    int const variable_id, int const component_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim)
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
            "there "
            "was no python script file specified.",
            bc_object);
    }

    auto* bc = scope[bc_object.c_str()]
                   .cast<PythonBoundaryConditionPythonSideInterface*>();

    if (variable_id >= static_cast<int>(dof_table.getNumberOfVariables()) ||
        component_id >= dof_table.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: ({:d}, "
            "{:d}), "
            "maximum values: ({:d}, {:d}).",
            variable_id, component_id, dof_table.getNumberOfVariables(),
            dof_table.getNumberOfVariableComponents(variable_id));
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
        PythonBoundaryConditionData{
            bc, dof_table, bulk_mesh_id,
            dof_table.getGlobalComponent(variable_id, component_id),
            boundary_mesh},
        integration_order, shapefunction_order, global_dim, flush_stdout);
}

}  // namespace ProcessLib
