/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/Algorithm.h"
#include "NumLib/Extrapolation/ExtrapolatableElementCollection.h"
#include "NumLib/Extrapolation/Extrapolator.h"
#include "ProcessLib/Utils/TransposeInPlace.h"

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
     * \remark The \c result_cache can be used to store the \c GlobalVector if
     * it is computed on-the-fly. Then a reference to the result cache can be
     * returned. Otherwise the \c Function must return a reference to a \c
     * GlobalVector that is stored somewhere else.
     */
    using Function = std::function<GlobalVector const&(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::unique_ptr<GlobalVector>& result_cache)>;

    template <typename F1, typename F2>
    SecondaryVariableFunctions(const unsigned num_components_, F1&& eval_field_,
                               F2&& eval_residuals_)
        : num_components(num_components_),
          eval_field(std::forward<F1>(eval_field_)),
          eval_residuals(std::forward<F2>(eval_residuals_))
    {
        // Used to detect nasty implicit conversions.
        static_assert(
            std::is_same_v<
                GlobalVector const&,
                typename std::invoke_result_t<
                    F1, double const, std::vector<GlobalVector*> const&,
                    std::vector<NumLib::LocalToGlobalIndexMap const*> const&,
                    std::unique_ptr<GlobalVector>&>>,
            "The function eval_field_ does not return a const reference"
            " to a GlobalVector");

        static_assert(
            std::is_same_v<
                GlobalVector const&,
                typename std::invoke_result_t<
                    F2, double const, std::vector<GlobalVector*> const&,
                    std::vector<NumLib::LocalToGlobalIndexMap const*> const&,
                    std::unique_ptr<GlobalVector>&>>,
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
            std::is_same_v<
                GlobalVector const&,
                typename std::invoke_result_t<
                    F1, double const, std::vector<GlobalVector*> const&,
                    std::vector<NumLib::LocalToGlobalIndexMap const*> const&,
                    std::unique_ptr<GlobalVector>&>>,
            "The function eval_field_ does not return a const reference"
            " to a GlobalVector");
    }

    const unsigned num_components;  //!< Number of components of the variable.

    //! Computes the value of the field at every node of the underlying mesh.
    Function const eval_field;

    //! If the secondary variable is extrapolated from integration points to
    //! mesh nodes, this function computes the extrapolation residual. For
    //! further information check the specific NumLib::Extrapolator
    //! documentation.
    Function const eval_residuals;
};

//! Stores information about a specific secondary variable
struct SecondaryVariable final
{
    std::string const name;  //!< Name of the variable; used, e.g., for output.

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
    SecondaryVariable const& get(std::string const& external_name) const;

    std::map<std::string, std::string>::const_iterator begin() const;
    std::map<std::string, std::string>::const_iterator end() const;

private:
    //! Maps external variable names to internal ones.
    //! The external variable names are used, e.g., for output.
    std::map<std::string, std::string> _map_external_to_internal;

    //! Collection of all configured secondary variables.
    //! Maps the internal variable name to the corresponding SecondaryVariable
    //! instance.
    std::map<std::string, SecondaryVariable> _configured_secondary_variables;
};

/*! Creates an object that computes a secondary variable via extrapolation of
 * integration point values.
 *
 * \param num_components The number of components of the secondary variable.
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
    auto const eval_field =
        [num_components, &extrapolator, &local_assemblers,
         integration_point_values_method](
            const double t,
            std::vector<GlobalVector*> const& x,
            std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
            std::unique_ptr<GlobalVector>& /*result_cache*/
            ) -> GlobalVector const&
    {
        auto const extrapolatables = NumLib::makeExtrapolatable(
            local_assemblers, integration_point_values_method);
        extrapolator.extrapolate(num_components, extrapolatables, t, x,
                                 dof_tables);
        return extrapolator.getNodalValues();
    };

    auto const eval_residuals =
        [num_components, &extrapolator, &local_assemblers,
         integration_point_values_method](
            const double t,
            std::vector<GlobalVector*> const& x,
            std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
            std::unique_ptr<GlobalVector>& /*result_cache*/
            ) -> GlobalVector const&
    {
        auto const extrapolatables = NumLib::makeExtrapolatable(
            local_assemblers, integration_point_values_method);
        extrapolator.calculateResiduals(num_components, extrapolatables, t, x,
                                        dof_tables);
        return extrapolator.getElementResiduals();
    };
    return {num_components, eval_field, eval_residuals};
}

//! A variant that takes an \c accessor function as a callback.
//!
//! The \c accessor must take a local assembler by const reference as its only
//! argument and must return a <tt>std::vector\<double\></tt> of integration
//! point data. The data returned by the accessor must be in "integration point
//! writer order", which is transposed compared to "extrapolator order".
template <typename LocalAssemblerCollection, typename IPDataAccessor>
SecondaryVariableFunctions makeExtrapolator2(
    const unsigned num_components,
    NumLib::Extrapolator& extrapolator,
    LocalAssemblerCollection const& local_assemblers,
    IPDataAccessor&& accessor)
{
    using LocalAssemblerInterface = std::remove_cvref_t<
        decltype(*std::declval<LocalAssemblerCollection>()[0])>;
    static_assert(std::is_invocable_r_v<std::vector<double>, IPDataAccessor,
                                        LocalAssemblerInterface const&>);

    if (num_components == 1)
    {
        auto method_wrapped =
            [accessor](
                LocalAssemblerInterface const& loc_asm, const double /*t*/,
                std::vector<GlobalVector*> const& /*x*/,
                std::vector<
                    NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
                std::vector<double>& cache) -> std::vector<double> const&
        {
            cache = accessor(loc_asm);
            return cache;
        };

        return makeExtrapolator(num_components, extrapolator, local_assemblers,
                                method_wrapped);
    }

    auto method_wrapped =
        [accessor, num_components](
            LocalAssemblerInterface const& loc_asm, const double /*t*/,
            std::vector<GlobalVector*> const& /*x*/,
            std::vector<
                NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
            std::vector<double>& cache) -> std::vector<double> const&
    {
        cache = accessor(loc_asm);
        transposeInPlace(cache, cache.size() / num_components);
        return cache;
    };

    return makeExtrapolator(num_components, extrapolator, local_assemblers,
                            method_wrapped);
}

}  // namespace ProcessLib
