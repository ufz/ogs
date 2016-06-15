/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SECONDARY_VARIABLE_H
#define PROCESSLIB_SECONDARY_VARIABLE_H

#include "BaseLib/ConfigTree.h"
#include "BaseLib/uniqueInsert.h"
#include "NumLib/Extrapolation/Extrapolator.h"

namespace NumLib { class LocalToGlobalIndexMap; }

namespace ProcessLib
{

//! Holder for function objects that compute secondary variables,
//! and (optionally) also the residuals (e.g., in case of extrapolation)
struct SecondaryVariableFunctions final
{
    /*! Type of functions used.
     *
     * \note The argument \c dof_table is the d.o.f. table of the process, i.e.
     * it possibly contains information about several process variables.
     *
     * \remark The \c result_cache can be used to store the \c GlobalVector if it
     * is computed on-the-fly. Then a reference to the result cache can be returned.
     * Otherwise the \c Function must return a reference to a \c GlobalVector that
     * is stored somewhere else.
     */
    using Function = std::function<GlobalVector const&(
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache)>;

    SecondaryVariableFunctions() = default;

    template<typename F1, typename F2>
    SecondaryVariableFunctions(F1&& eval_field_, F2&& eval_residuals_)
        : eval_field(std::forward<F1>(eval_field_))
        , eval_residuals(std::forward<F2>(eval_residuals_))
    {
        // Used to detect nasty implicit conversions.
        static_assert(std::is_same<GlobalVector const&,
            typename std::result_of<F1(
                GlobalVector const&, NumLib::LocalToGlobalIndexMap const&,
                std::unique_ptr<GlobalVector>&
                )>::type>::value,
            "The function eval_field_ does not return a const reference"
            " to a GlobalVector");

        static_assert(std::is_same<GlobalVector const&,
            typename std::result_of<F2(
                GlobalVector const&, NumLib::LocalToGlobalIndexMap const&,
                std::unique_ptr<GlobalVector>&
            )>::type>::value,
            "The function eval_residuals_ does not return a const reference"
            " to a GlobalVector");
    }

    template<typename F1>
    SecondaryVariableFunctions(F1&& eval_field_, std::nullptr_t)
        : eval_field(std::forward<F1>(eval_field_))
    {
        // Used to detect nasty implicit conversions.
        static_assert(std::is_same<GlobalVector const&,
            typename std::result_of<F1(
                GlobalVector const&, NumLib::LocalToGlobalIndexMap const&,
                std::unique_ptr<GlobalVector>&
                )>::type>::value,
            "The function eval_field_ does not return a const reference"
            " to a GlobalVector");
    }

    Function eval_field;
    Function eval_residuals;
};

//! Stores information about a specific secondary variable
struct SecondaryVariable final
{
    std::string const name;      //!< Name of the variable; used, e.g., for output.
    const unsigned n_components; //!< Number of components of the variable.

    //! Functions used for computing the secondary variable.
    SecondaryVariableFunctions fcts;
};

//! Handles configuration of several secondary variables from the project file.
class SecondaryVariableCollection final
{
public:
    /*! Constructs new instance.
     *
     * \param config    Configuration settings.
     * \param tag_names Possible tag names that contain information about specific
     *                  secondary variables.
     *
     * In this method only the mapping from tag names to variables is set up.
     * Any further information has to be passed via addSecondaryVariable().
     */
    SecondaryVariableCollection(
            boost::optional<BaseLib::ConfigTree> const& config,
            std::initializer_list<std::string> tag_names)
    {
        if (!config) return;

        // read which variables are defined in the config
        for (auto const& tag_name : tag_names) {
            if (!_all_secondary_variables.insert(tag_name).second) {
                OGS_FATAL("Tag name <%s> has been specified twice as a secondary variable.");
            }

            //! \ogs_file_special
            if (auto var_name = config->getConfigParameterOptional<std::string>(tag_name))
            {
                // TODO check primary vars, too
                BaseLib::insertIfKeyValueUniqueElseError(
                            _map_tagname_to_varname, tag_name, *var_name,
                            "Secondary variable names must be unique.");
            }
        }
    }

    /*! Tells if a secondary variable with the specified name has been set up.
     *
     * \note \c variable_name is not the tag name in the project file!
     */
    bool variableExists(std::string const& variable_name) const
    {
        auto pred = [&variable_name](std::pair<std::string, std::string> const& p) {
            return p.second == variable_name;
        };

        // check if out_var is a  secondary variable
        auto const& var = std::find_if(
            _map_tagname_to_varname.cbegin(), _map_tagname_to_varname.cend(), pred);

        return var != _map_tagname_to_varname.cend();
    }

    /*! Set up a secondary variable.
     *
     * \param tag_name the tag in the project file associated with this
     * secondary variable.
     * \param num_components the variable's number of components.
     * \param fcts functions that compute the variable.
     *
     * \note
     * Only variables requested by the user in the project file will be
     * configured.
     * All other variables are silently ignored.
     */
    void addSecondaryVariable(std::string const& tag_name,
                              const unsigned num_components,
                              SecondaryVariableFunctions&& fcts)
    {
        auto it = _map_tagname_to_varname.find(tag_name);

        // get user-supplied var_name for the given tag_name
        if (it != _map_tagname_to_varname.end())
        {
            auto const& var_name = it->first;

            if (!_configured_secondary_variables
                     .emplace(std::make_pair(
                         var_name,
                         SecondaryVariable{
                             var_name, num_components, std::move(fcts)}))
                     .second)
            {
                OGS_FATAL("The secondary variable with name `%s' has already been "
                          "set up.",
                          var_name.c_str());
            }
        }
        else if (_all_secondary_variables.find(tag_name) ==
                 _all_secondary_variables.end())
        {
            OGS_FATAL("The tag <%s> has not been registered to mark a secondary "
                      "variable.",
                      tag_name.c_str());
        }
    }

    //! Returns an iterator to the first secondary variable.
    typename std::map<std::string,
                      SecondaryVariable>::const_iterator
    begin() const
    {
        return _configured_secondary_variables.begin();
    }

    //! Returns an iterator past the last secondary variable.
    typename std::map<std::string,
                      SecondaryVariable>::const_iterator
    end() const
    {
        return _configured_secondary_variables.end();
    }

private:
    //! Maps project file tag names to secondary variable names.
    std::map<std::string, std::string> _map_tagname_to_varname;

    //! Collection of all configured secondary variables.
    //! Maps the variable name to the corresponding SecondaryVariable.
    std::map<std::string, SecondaryVariable> _configured_secondary_variables;

    //! Set of all tags available as a secondary variable.
    std::set<std::string> _all_secondary_variables;
};


//! Creates an object that computes a secondary variable via extrapolation
//! of integration point values.
template<typename PropertyEnum, typename LocalAssembler>
SecondaryVariableFunctions
makeExtrapolator(PropertyEnum const property,
                 NumLib::Extrapolator<PropertyEnum, LocalAssembler>&
                 extrapolator,
                 typename NumLib::Extrapolator<PropertyEnum,
                     LocalAssembler>::LocalAssemblers const& local_assemblers)
{
    static_assert(std::is_base_of<
         NumLib::Extrapolatable<PropertyEnum>, LocalAssembler>::value,
        "The passed local assembler type (i.e. the local assembler interface) must"
        " derive from NumLib::Extrapolatable<>.");

    auto const eval_field = [property, &extrapolator, &local_assemblers](
            GlobalVector const& /*x*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
            std::unique_ptr<GlobalVector>& /*result_cache*/
            ) -> GlobalVector const&
    {
        extrapolator.extrapolate(local_assemblers, property);
        return extrapolator.getNodalValues();
    };

    auto const eval_residuals = [property, &extrapolator, &local_assemblers](
            GlobalVector const& /*x*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
            std::unique_ptr<GlobalVector>& /*result_cache*/
            ) -> GlobalVector const&
    {
        extrapolator.calculateResiduals(local_assemblers, property);
        return extrapolator.getElementResiduals();
    };
    return { eval_field, eval_residuals };
}

} // namespace ProcessLib

#endif // PROCESSLIB_SECONDARY_VARIABLE_H
