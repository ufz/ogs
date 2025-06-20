/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstraintDirichletBoundaryCondition.h"

#include <algorithm>
#include <vector>

#include "BaseLib/Logging.h"
#include "MeshLib/MeshSearch/NodeSearch.h"  // for getUniqueNodes
#include "MeshLib/Node.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
ConstraintDirichletBoundaryCondition::ConstraintDirichletBoundaryCondition(
    ParameterLib::Parameter<double> const& parameter,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id, MeshLib::Mesh const& bc_mesh,
    unsigned const integration_order, MeshLib::Mesh const& bulk_mesh,
    double const constraint_threshold, bool const lower,
    std::function<Eigen::Vector3d(std::size_t const, MathLib::Point3d const&,
                                  double const,
                                  std::vector<GlobalVector*> const&)>
        getFlux)
    : _parameter(parameter),
      _variable_id(variable_id),
      _component_id(component_id),
      _bc_mesh(bc_mesh),
      _integration_order(integration_order),
      _constraint_threshold(constraint_threshold),
      _lower(lower),
      _bulk_mesh(bulk_mesh),
      _getFlux(getFlux)
{
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

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    DBUG(
        "Found {:d} nodes for constraint Dirichlet BCs for the variable {:d} "
        "and component {:d}",
        bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset{_bc_mesh, bc_nodes};

    // Create local DOF table from intersected mesh subsets for the given
    // variable and component ids.
    _dof_table_boundary = dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset));

    auto const& bc_elements(_bc_mesh.getElements());
    _local_assemblers.resize(bc_elements.size());
    _flux_values.resize(bc_elements.size());
    // create _bulk_ids vector
    auto const* bulk_element_ids = MeshLib::bulkElementIDs(_bc_mesh);
    auto const* bulk_node_ids = MeshLib::bulkNodeIDs(_bc_mesh);
    auto const& bulk_nodes = bulk_mesh.getNodes();

    auto get_bulk_element_face_id =
        [&](auto const bulk_element_id, MeshLib::Element const* bc_elem)
    {
        auto const* bulk_elem = _bulk_mesh.getElement(bulk_element_id);
        std::array<MeshLib::Node const*, 3> nodes{
            {bulk_nodes[(*bulk_node_ids)[bc_elem->getNode(0)->getID()]],
             bulk_nodes[(*bulk_node_ids)[bc_elem->getNode(1)->getID()]],
             bulk_nodes[(*bulk_node_ids)[bc_elem->getNode(2)->getID()]]}};
        return bulk_elem->identifyFace(nodes.data());
    };

    _bulk_ids.reserve(bc_elements.size());
    std::transform(
        begin(bc_elements), end(bc_elements), std::back_inserter(_bulk_ids),
        [&](auto const* bc_element)
        {
            auto const bulk_element_id =
                (*bulk_element_ids)[bc_element->getID()];
            return std::make_pair(
                bulk_element_id,
                get_bulk_element_face_id(bulk_element_id, bc_element));
        });

    const int shape_function_order = 1;

    BoundaryConditionAndSourceTerm::createLocalAssemblers<
        ConstraintDirichletBoundaryConditionLocalAssembler>(
        _bulk_mesh.getDimension(), _bc_mesh.getElements(), *_dof_table_boundary,
        shape_function_order, _local_assemblers,
        NumLib::IntegrationOrder{_integration_order},
        _bc_mesh.isAxiallySymmetric(), bulk_mesh, _bulk_ids);
}

void ConstraintDirichletBoundaryCondition::preTimestep(
    double const t, std::vector<GlobalVector*> const& x,
    int const /*process_id*/)
{
    DBUG(
        "ConstraintDirichletBoundaryCondition::preTimestep: computing flux "
        "constraints");
    for (auto const id : _bc_mesh.getElements() | MeshLib::views::ids)
    {
        _flux_values[id] = _local_assemblers[id]->integrate(
            x, t,
            [this](std::size_t const element_id, MathLib::Point3d const& pnt,
                   double const t, std::vector<GlobalVector*> const& x)
            { return _getFlux(element_id, pnt, t, x); });
    }
}

void ConstraintDirichletBoundaryCondition::getEssentialBCValues(
    const double t, const GlobalVector& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    bc_values.ids.clear();
    bc_values.values.clear();

    std::vector<std::pair<GlobalIndexType, double>> tmp_bc_values;

    auto isFlux = [&](const std::size_t element_id)
    {
        return _lower ? _flux_values[element_id] < _constraint_threshold
                      : _flux_values[element_id] > _constraint_threshold;
    };

    for (auto const* boundary_element : _bc_mesh.getElements())
    {
        // check if the boundary element is active
        if (isFlux(boundary_element->getID()))
        {
            continue;
        }

        // loop over the boundary element nodes and set the dirichlet value
        unsigned const number_nodes = boundary_element->getNumberOfNodes();
        for (unsigned i = 0; i < number_nodes; ++i)
        {
            auto const& node = *boundary_element->getNode(i);
            auto const id = node.getID();
            ParameterLib::SpatialPosition const pos{
                id, boundary_element->getID(), node};

            MeshLib::Location l(_bc_mesh.getID(), MeshLib::MeshItemType::Node,
                                id);
            // TODO: that might be slow, but only done once
            const auto g_idx = _dof_table_boundary->getGlobalIndex(
                l, _variable_id, _component_id);
            if (g_idx == NumLib::MeshComponentMap::nop)
            {
                continue;
            }
            // For the DDC approach (e.g. with PETSc option), the negative
            // index of g_idx means that the entry by that index is a ghost one,
            // which should be dropped. Especially for PETSc routines
            // MatZeroRows and MatZeroRowsColumns, which are called to apply the
            // Dirichlet BC, the negative index is not accepted like other
            // matrix or vector PETSc routines. Therefore, the following
            // if-condition is applied.
            if (g_idx >= 0)
            {
                tmp_bc_values.emplace_back(g_idx, _parameter(t, pos).front());
            }
        }
    }

    if (tmp_bc_values.empty())
    {
        DBUG("The domain on which the boundary is defined is empty");
        return;
    }

    // store unique pairs of node id and value with respect to the node id in
    // the bc_values vector. The values related to the same node id are
    // averaged.
    // first: sort the (node id, value) pairs according to the node id
    std::sort(tmp_bc_values.begin(), tmp_bc_values.end(),
              [](std::pair<GlobalIndexType, double> const& a,
                 std::pair<GlobalIndexType, double> const& b)
              { return a.first < b.first; });
    // second: average the values over equal node id ranges
    unsigned cnt = 1;
    GlobalIndexType current_id = tmp_bc_values.begin()->first;
    double sum = tmp_bc_values.begin()->second;
    for (auto const& tmp_bc_value : tmp_bc_values)
    {
        if (tmp_bc_value.first == current_id)
        {
            cnt++;
            sum += tmp_bc_value.second;
        }
        else
        {
            bc_values.ids.emplace_back(current_id);
            bc_values.values.emplace_back(sum / cnt);
            cnt = 1;
            sum = tmp_bc_value.second;
            current_id = tmp_bc_value.first;
        }
    }
    bc_values.ids.emplace_back(current_id);
    bc_values.values.emplace_back(sum / cnt);

    DBUG("Found {:d} constraint dirichlet boundary condition values.",
         bc_values.ids.size());
}

std::unique_ptr<ConstraintDirichletBoundaryCondition>
createConstraintDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    unsigned const integration_order, int const component_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    Process const& constraining_process)
{
    DBUG("Constructing ConstraintDirichletBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "ConstraintDirichlet");

    auto const constraint_type =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__ConstraintDirichlet__constraint_type}
        config.getConfigParameter<std::string>("constraint_type");
    if (constraint_type != "Flux")
    {
        OGS_FATAL("The constraint type is '{:s}', but has to be 'Flux'.",
                  constraint_type);
    }

    // Todo (TF) Open question: How to specify which getFlux function should be
    // used for the constraint calculation?
    auto const constraining_process_variable =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__ConstraintDirichlet__constraining_process_variable}
        config.getConfigParameter<std::string>("constraining_process_variable");

    if (!constraining_process.isMonolithicSchemeUsed())
    {
        OGS_FATAL(
            "The constraint dirichlet boundary condition is implemented only "
            "for monolithic implemented processes.");
    }
    const int process_id = 0;
    auto process_variables =
        constraining_process.getProcessVariables(process_id);
    auto constraining_pv =
        std::find_if(process_variables.cbegin(), process_variables.cend(),
                     [&constraining_process_variable](ProcessVariable const& pv)
                     { return pv.getName() == constraining_process_variable; });
    if (constraining_pv == std::end(process_variables))
    {
        auto const& constraining_process_variable_name =
            process_variables[variable_id].get().getName();
        OGS_FATAL(
            "<constraining_process_variable> in process variable name '{:s}' "
            "at geometry 'TODO' : The constraining process variable is set as "
            "'{:s}', but this is not specified in the project file.",
            constraining_process_variable_name,
            constraining_process_variable);
    }

    auto const constraint_threshold =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__ConstraintDirichlet__constraint_threshold}
        config.getConfigParameter<double>("constraint_threshold");

    auto const constraint_direction_string =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__ConstraintDirichlet__constraint_direction}
        config.getConfigParameter<std::string>("constraint_direction");
    if (constraint_direction_string != "greater" &&
        constraint_direction_string != "lower")
    {
        OGS_FATAL(
            "The constraint direction is '{:s}', but has to be either "
            "'greater' or 'lower'.",
            constraint_direction_string);
    }
    bool const lower = constraint_direction_string == "lower";

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__ConstraintDirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter {:s}", param_name);

    auto& param = ParameterLib::findParameter<double>(param_name, parameters, 1,
                                                      &bc_mesh);

    // In case of partitioned mesh the boundary could be empty, i.e. there is no
    // boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createBoundaryCondition().
    if (bc_mesh.getDimension() == 0 && bc_mesh.getNumberOfNodes() == 0 &&
        bc_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<ConstraintDirichletBoundaryCondition>(
        param, dof_table_bulk, variable_id, component_id, bc_mesh,
        integration_order, constraining_process.getMesh(), constraint_threshold,
        lower,
        [&constraining_process](std::size_t const element_id,
                                MathLib::Point3d const& pnt, double const t,
                                std::vector<GlobalVector*> const& x)
        { return constraining_process.getFlux(element_id, pnt, t, x); });
}

}  // namespace ProcessLib
