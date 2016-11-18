/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "BaseLib/uniqueInsert.h"
#include "NumLib/Extrapolation/Extrapolator.h"
#include "NumLib/Extrapolation/ExtrapolatableElementCollection.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}

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
        const double t,
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache)>;

    SecondaryVariableFunctions() = default;

    template <typename F1, typename F2>
    SecondaryVariableFunctions(const unsigned num_components_, F1&& eval_field_,
                               F2&& eval_residuals_)
        : num_components(num_components_),
          eval_field(std::forward<F1>(eval_field_)),
          eval_residuals(std::forward<F2>(eval_residuals_))
    {
        // Used to detect nasty implicit conversions.
        static_assert(
            std::is_same<GlobalVector const&,
                         typename std::result_of<F1(
                             double const, GlobalVector const&,
                             NumLib::LocalToGlobalIndexMap const&,
                             std::unique_ptr<GlobalVector>&)>::type>::value,
            "The function eval_field_ does not return a const reference"
            " to a GlobalVector");

        static_assert(
            std::is_same<GlobalVector const&,
                         typename std::result_of<F2(
                             double const, GlobalVector const&,
                             NumLib::LocalToGlobalIndexMap const&,
                             std::unique_ptr<GlobalVector>&)>::type>::value,
            "The function eval_residuals_ does not return a const reference"
            " to a GlobalVector");
    }

    template <typename F1>
    SecondaryVariableFunctions(const unsigned num_components_, F1&& eval_field_,
                               std::nullptr_t)
        : num_components(num_components_),
          eval_field(std::forward<F1>(eval_field_))
    {
        // Used to detect nasty implicit conversions.
        static_assert(
            std::is_same<GlobalVector const&,
                         typename std::result_of<F1(
                             double const, GlobalVector const&,
                             NumLib::LocalToGlobalIndexMap const&,
                             std::unique_ptr<GlobalVector>&)>::type>::value,
            "The function eval_field_ does not return a const reference"
            " to a GlobalVector");
    }

    const unsigned num_components;  //!< Number of components of the variable.
    Function const eval_field;
    Function const eval_residuals;
};

//! Stores information about a specific secondary variable
struct SecondaryVariable final
{
    std::string const name;      //!< Name of the variable; used, e.g., for output.

    //! Functions used for computing the secondary variable.
    SecondaryVariableFunctions fcts;
};

//! Handles configuration of several secondary variables from the project file.
class SecondaryVariableCollection final
{
public:
    //! Register a variable with the given internal and external names.
    void addNameMapping(std::string const& internal_name,
                        std::string const& external_name);

    /*! Set up a secondary variable.
     *
     * \param internal_name the tag in the project file associated with this
     * secondary variable.
     * \param fcts functions that compute the variable.
     *
     * \note
     * Only variables requested by the user in the project file will be
     * configured.
     * All other variables are silently ignored.
     */
    void addSecondaryVariable(std::string const& internal_name,
                              SecondaryVariableFunctions&& fcts);

    //! Returns the secondary variable with the given external name.
    SecondaryVariable const& get(std::string const& external_name);

private:
    //! Maps external variable names to internal ones.
    //! The external variable names are used, e.g., for output.
    std::map<std::string, std::string> _map_external_to_internal;

    //! Collection of all configured secondary variables.
    //! Maps the internal variable name to the corresponding SecondaryVariable
    //! instance.
    std::map<std::string, SecondaryVariable> _configured_secondary_variables;

    //! Set of all internal variable names known to this instance.
    std::set<std::string> _all_secondary_variables;
};

/*! Creates an object that computes a secondary variable via extrapolation of
 * integration point values.
 *
 * \param extrapolator The extrapolator used for extrapolation.
 * \param local_assemblers The collection of local assemblers whose integration
 * point values will be extrapolated.
 * \param integration_point_values_method The member function of the local
 * assembler returning/computing the integration point values of the specific
 * property being extrapolated.
 */
template <typename LocalAssemblerCollection>
SecondaryVariableFunctions makeExtrapolator(
    const unsigned num_components,
    NumLib::Extrapolator& extrapolator,
    LocalAssemblerCollection const& local_assemblers,
    typename NumLib::ExtrapolatableLocalAssemblerCollection<
        LocalAssemblerCollection>::IntegrationPointValuesMethod
        integration_point_values_method)
{
    auto const eval_field = [num_components, &extrapolator, &local_assemblers,
                             integration_point_values_method](
        const double t,
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector> & /*result_cache*/
        ) -> GlobalVector const& {
        auto const extrapolatables = NumLib::makeExtrapolatable(
            local_assemblers, integration_point_values_method);
        extrapolator.extrapolate(num_components, extrapolatables, t, x,
                                 dof_table);
        return extrapolator.getNodalValues();
    };

    auto const eval_residuals = [num_components, &extrapolator,
                                 &local_assemblers,
                                 integration_point_values_method](
        const double t,
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector> & /*result_cache*/
        ) -> GlobalVector const& {
        auto const extrapolatables = NumLib::makeExtrapolatable(
            local_assemblers, integration_point_values_method);
        extrapolator.calculateResiduals(num_components, extrapolatables, t, x,
                                        dof_table);
        return extrapolator.getElementResiduals();
    };
    return {num_components, eval_field, eval_residuals};
}

}  // namespace ProcessLib
