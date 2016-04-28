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
#include "NumLib/Extrapolation/Extrapolator.h"
#include "ProcessVariable.h"

namespace ProcessLib
{

template<typename GlobalVector>
struct SecondaryVariableFunctions
{
    using Fct = std::function<GlobalVector(
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
            ) -> GlobalVector
    {
        extrapolator.extrapolate(local_assemblers, property);
        return extrapolator.getNodalValues();
    };

    auto const eval_residuals = [property, &extrapolator, &local_assemblers](
            GlobalVector const& /*x*/,
            AssemblerLib::LocalToGlobalIndexMap const& /*dof_table*/
            ) -> GlobalVector
    {
        extrapolator.calculateResiduals(local_assemblers, property);
        return extrapolator.getElementResiduals();
    };
    return { eval_field, eval_residuals };
}


template<typename GlobalVector>
struct SecondaryVariable
{
    std::string const name;
    const unsigned n_components;
    SecondaryVariableFunctions<GlobalVector> fcts;
};

template <typename GlobalVector>
struct ProcessOutput
{
    ProcessOutput() = default;
    ProcessOutput(ProcessOutput<GlobalVector>&&) = default;

    //! There is no need to copy instances of this struct.
    ProcessOutput(ProcessOutput<GlobalVector> const&) = delete;

    void addSecondaryVariable(
            BaseLib::ConfigTree const& config,
            std::string const& var_tag, const unsigned num_components,
            SecondaryVariableFunctions<GlobalVector>&& fcts)
    {
        if (auto var_name = config.getConfParamOptional<std::string>(var_tag))
        {
            secondary_variables.push_back(
                {*var_name, num_components, std::move(fcts)});
        }
    }

    template<typename ProcVarRef>
    void setOutputVariables(BaseLib::ConfigTree const& output_config,
                            std::vector<ProcVarRef> const& process_variables)
    {
        auto out_vars = output_config.getConfSubtreeOptional("variables");
        if (!out_vars) return;

        for (auto out_var : out_vars->getConfParamList<std::string>("variable"))
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

            if (pcs_var == process_variables.cend())
            {
                auto pred2 = [&out_var](SecondaryVariable<GlobalVector> const& p) {
                    return p.name == out_var;
                };

                // check if out_var is a  secondary variable
                auto const& pcs_var2 = std::find_if(
                    secondary_variables.cbegin(), secondary_variables.cend(), pred2);

                if (pcs_var2 == secondary_variables.cend())
                {
                    ERR("Output variable `%s' is neither a process variable nor a"
                        " secondary variable", out_var.c_str());
                    std::abort();
                }
            }

            DBUG("adding output variable `%s'", out_var.c_str());
            output_variables.insert(out_var);
        }

        if (auto out_resid = output_config.getConfParamOptional<bool>(
                "output_extrapolation_residuals")) {
            output_residuals = *out_resid;
        }
    }

    std::vector<SecondaryVariable<GlobalVector>> secondary_variables;
    std::set<std::string> output_variables;

    bool output_residuals = false;
    //! Output global matrix/rhs after first iteration.
    bool output_global_matrix = false;
    bool output_iteration_results = false;
};

} // ProcessLib


#endif // PROCESSLIB_PROCESSOUTPUT_H
