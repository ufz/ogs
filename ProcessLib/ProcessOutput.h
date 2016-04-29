/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_PROCESSOUTPUT_H
#define PROCESSLIB_PROCESSOUTPUT_H

#include "BaseLib/ConfigTree.h"
#include "BaseLib/uniqueInsert.h"
#include "NumLib/Extrapolation/Extrapolator.h"
#include "ProcessVariable.h"

namespace ProcessLib
{

template<typename GlobalVector>
struct SecondaryVariableFunctions final
{
    using Fct = std::function<GlobalVector const&(
        GlobalVector const& x,
        AssemblerLib::LocalToGlobalIndexMap const& dof_table)>;

    Fct eval_field;
    Fct eval_residuals;
};

template<typename GlobalVector, typename PropertyEnum, typename LocalAssembler>
SecondaryVariableFunctions<GlobalVector>
makeExtrapolator(PropertyEnum const property,
                 NumLib::Extrapolator<GlobalVector, PropertyEnum, LocalAssembler>&
                 extrapolator,
                 typename NumLib::Extrapolator<GlobalVector, PropertyEnum,
                     LocalAssembler>::LocalAssemblers const& local_assemblers)
{
    static_assert(std::is_base_of<
         NumLib::Extrapolatable<GlobalVector, PropertyEnum>, LocalAssembler>::value,
        "The passed local assembler type (i.e. the local assembler interface) must"
        " derive from NumLib::Extrapolatable<>!");

    auto const eval_field = [property, &extrapolator, &local_assemblers](
            GlobalVector const& /*x*/,
            AssemblerLib::LocalToGlobalIndexMap const& /*dof_table*/
            ) -> GlobalVector const&
    {
        extrapolator.extrapolate(local_assemblers, property);
        return extrapolator.getNodalValues();
    };

    auto const eval_residuals = [property, &extrapolator, &local_assemblers](
            GlobalVector const& /*x*/,
            AssemblerLib::LocalToGlobalIndexMap const& /*dof_table*/
            ) -> GlobalVector const&
    {
        extrapolator.calculateResiduals(local_assemblers, property);
        return extrapolator.getElementResiduals();
    };
    return { eval_field, eval_residuals };
}


template<typename GlobalVector>
struct SecondaryVariable final
{
    std::string const name;
    const unsigned n_components;
    SecondaryVariableFunctions<GlobalVector> fcts;
};


template<typename GlobalVector>
class SecondaryVariableCollection final
{
public:
    SecondaryVariableCollection(
            boost::optional<BaseLib::ConfigTree> const& config,
            std::initializer_list<std::string> tag_names)
    {
        if (!config) return;

        // read which variables are defined in the config
        for (auto const& tag_name : tag_names) {
            if (auto var_name = config->getConfParamOptional<std::string>(tag_name))
            {
                // TODO check primary vars, too
                BaseLib::insertIfKeyValueUniqueElseError(
                            _map_tagname_to_varname, tag_name, *var_name,
                            "Secondary variable names must be unique.");
            }
        }
    }

    bool variableExists(std::string const& variable_name) const
    {
        auto pred = [&variable_name](SecondaryVariable<GlobalVector> const& p) {
            return p.name == variable_name;
        };

        // check if out_var is a  secondary variable
        auto const& var = std::find_if(
            _secondary_variables.cbegin(), _secondary_variables.cend(), pred);

        return var != _secondary_variables.cend();
    }

    void addSecondaryVariable(
            std::string const& tag_name, const unsigned num_components,
            SecondaryVariableFunctions<GlobalVector>&& fcts)
    {
        auto it = _map_tagname_to_varname.find(tag_name);

        // get user-supplied var_name for the given tag_name
        if (it != _map_tagname_to_varname.end()) {
            auto const& var_name = it->first;
            // TODO make sure the same variable is not pushed twice
            _secondary_variables.push_back(
                {var_name, num_components, std::move(fcts)});
        }
    }

    typename std::vector<SecondaryVariable<GlobalVector>>::const_iterator
    begin() const
    {
        return _secondary_variables.begin();
    }

    typename std::vector<SecondaryVariable<GlobalVector>>::const_iterator
    end() const
    {
        return _secondary_variables.end();
    }

private:
    std::map<std::string, std::string> _map_tagname_to_varname;
    std::vector<SecondaryVariable<GlobalVector>> _secondary_variables;
};


template <typename GlobalVector>
struct ProcessOutput final
{
    ProcessOutput(BaseLib::ConfigTree const& output_config,
                  std::vector<std::reference_wrapper<ProcessVariable>> const&
                  process_variables,
                  SecondaryVariableCollection<GlobalVector> const& secondary_variables)
    {
        auto const out_vars = output_config.getConfSubtree("variables");

        for (auto out_var : out_vars.getConfParamList<std::string>("variable"))
        {
            if (output_variables.find(out_var) != output_variables.cend())
            {
                ERR("output variable `%s' specified more than once.", out_var.c_str());
                std::abort();
            }

            auto pred = [&out_var](ProcessVariable const& pv) {
                return pv.getName() == out_var;
            };

            // check if out_var is a process variable
            auto const& pcs_var = std::find_if(
                process_variables.cbegin(), process_variables.cend(), pred);

            if (pcs_var == process_variables.cend()
                && !secondary_variables.variableExists(out_var))
            {
                ERR("Output variable `%s' is neither a process variable nor a"
                    " secondary variable", out_var.c_str());
                std::abort();
            }

            DBUG("adding output variable `%s'", out_var.c_str());
            output_variables.insert(out_var);
        }

        if (auto out_resid = output_config.getConfParamOptional<bool>(
                "output_extrapolation_residuals")) {
            output_residuals = *out_resid;
        }
    }

    std::set<std::string> output_variables;

    bool output_residuals = false;
    //! Output global matrix/rhs after first iteration.
    bool output_global_matrix = false;
    bool output_iteration_results = false;
};


template <typename GlobalVector>
void doProcessOutput(
        std::string const& file_name,
        GlobalVector const& x,
        MeshLib::Mesh& mesh,
        AssemblerLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
        SecondaryVariableCollection<GlobalVector> secondary_variables,
        ProcessOutput<GlobalVector> const& process_output)
{
    DBUG("Process output.");

    // Copy result
#ifdef USE_PETSC
    // TODO It is also possible directly to copy the data for single process
    // variable to a mesh property. It needs a vector of global indices and
    // some PETSc magic to do so.
    std::vector<double> x_copy(x.getLocalSize() + x.getGhostSize());
#else
    std::vector<double> x_copy(x.size());
#endif
    x.copyValues(x_copy);

    auto const& output_variables = process_output.output_variables;

    std::size_t const n_nodes = mesh.getNNodes();
    int global_component_offset = 0;
    int global_component_offset_next = 0;

    // primary variables
    for (ProcessVariable& pv : process_variables)
    {
        int const n_components = pv.getNumberOfComponents();
        global_component_offset = global_component_offset_next;
        global_component_offset_next += n_components;

        if (output_variables.find(pv.getName()) == output_variables.cend())
            continue;

        DBUG("  process variable %s", pv.getName().c_str());

        auto& output_data = pv.getOrCreateMeshProperty();

        for (std::size_t node_id = 0; node_id < n_nodes; ++node_id)
        {
            MeshLib::Location const l(mesh.getID(),
                                      MeshLib::MeshItemType::Node, node_id);
            // TODO extend component ids to multiple process variables.
            for (int component_id = 0; component_id < n_components;
                 ++component_id)
            {
                auto const global_component_id = global_component_offset + component_id;
                auto const index =
                        dof_table.getLocalIndex(
                            l, global_component_id, x.getRangeBegin(),
                            x.getRangeEnd());

                output_data[node_id * n_components + component_id] =
                        x_copy[index];
            }
        }
    }

    // the following section is for the output of secondary variables

    auto count_mesh_items = [](
        MeshLib::Mesh const& mesh, MeshLib::MeshItemType type) -> std::size_t
    {
        switch (type) {
        case MeshLib::MeshItemType::Cell: return mesh.getNElements();
        case MeshLib::MeshItemType::Node: return mesh.getNNodes();
        default: break; // avoid compiler warning
        }
        return 0;
    };

    auto get_or_create_mesh_property = [&mesh, &count_mesh_items](
        std::string const& property_name, MeshLib::MeshItemType type)
    {
        // Get or create a property vector for results.
        boost::optional<MeshLib::PropertyVector<double>&> result;

        auto const N = count_mesh_items(mesh, type);

        if (mesh.getProperties().hasPropertyVector(property_name))
        {
            result = mesh.getProperties().template
                     getPropertyVector<double>(property_name);
        }
        else
        {
            result = mesh.getProperties().template
                     createNewPropertyVector<double>(property_name, type);
            result->resize(N);
        }
        assert(result && result->size() == N);

        return result;
    };

    auto add_secondary_var = [&](SecondaryVariable<GlobalVector> const& var)
    {
        assert(var.n_components == 1); // TODO implement other cases

        {
            DBUG("  secondary variable %s", var.name.c_str());

            auto result = get_or_create_mesh_property(
                              var.name, MeshLib::MeshItemType::Node);
            assert(result->size() == mesh.getNNodes());

            auto const& nodal_values =
                    var.fcts.eval_field(x, dof_table);

            // Copy result
            for (std::size_t i = 0; i < mesh.getNNodes(); ++i)
            {
                assert(!std::isnan(nodal_values[i]));
                (*result)[i] = nodal_values[i];
            }
        }

        if (process_output.output_residuals && var.fcts.eval_residuals)
        {
            DBUG("  process var %s residual", var.name.c_str());
            auto const& property_name_res = var.name + "_residual";

            auto result = get_or_create_mesh_property(
                              property_name_res, MeshLib::MeshItemType::Cell);
            assert(result->size() == mesh.getNElements());

            auto const& residuals =
                    var.fcts.eval_residuals(x, dof_table);

            // Copy result
            for (std::size_t i = 0; i < mesh.getNElements(); ++i)
            {
                assert(!std::isnan(residuals[i]));
                (*result)[i] = residuals[i];
            }
        }
    };

    for (auto const& var : secondary_variables)
    {
        if (output_variables.find(var.name) != output_variables.cend())
            add_secondary_var(var);
    }

    // secondary variables output end

    // Write output file
    DBUG("Writing output to \'%s\'.", file_name.c_str());
    FileIO::VtuInterface vtu_interface(&mesh, vtkXMLWriter::Binary, true);
    vtu_interface.writeToFile(file_name);
}

} // ProcessLib


#endif // PROCESSLIB_PROCESSOUTPUT_H
