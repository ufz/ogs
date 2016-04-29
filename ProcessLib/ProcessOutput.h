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

} // ProcessLib


#endif // PROCESSLIB_PROCESSOUTPUT_H
