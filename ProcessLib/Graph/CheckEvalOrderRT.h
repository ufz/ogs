// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <boost/mp11.hpp>
#include <typeindex>
#include <unordered_set>

#include "Apply.h"
#include "BaseLib/DemangleTypeInfo.h"
#include "BaseLib/Logging.h"

namespace ProcessLib::Graph
{
namespace detail
{
template <typename T>
struct IsInputArgument
    : boost::mp11::mp_bool<std::is_lvalue_reference_v<T> &&
                           std::is_const_v<std::remove_reference_t<T>>>
{
    static_assert(std::is_lvalue_reference_v<T>,
                  "The current implementation only deals with l-value "
                  "references as function arguments. If you want to extend it, "
                  "test thoroughly in order to not introduce bugs.");
};

template <typename T>
struct IsOutputArgument : boost::mp11::mp_bool<!IsInputArgument<T>::value>
{
    static_assert(std::is_lvalue_reference_v<T>,
                  "The current implementation only deals with l-value "
                  "references as function arguments. If you want to extend it, "
                  "test thoroughly in order to not introduce bugs.");
};

template <typename Model>
bool isEvalOrderCorrectRT(std::unordered_set<std::type_index>& computed_data)
{
    using namespace boost::mp11;

    using ModelArgs =
        typename GetFunctionArgumentTypes<decltype(&Model::eval)>::type;
    using ModelInputs = mp_filter<IsInputArgument, ModelArgs>;
    using ModelOutputs = mp_filter<IsOutputArgument, ModelArgs>;

    using ModelInputsWrapped = mp_transform<mp_identity, ModelInputs>;

    // Check that all inputs have already been computed before.
    bool all_inputs_computed = true;
    mp_for_each<ModelInputsWrapped>(
        [&computed_data,
         &all_inputs_computed]<typename Input>(mp_identity<Input>)
        {
            if (!computed_data.contains(std::type_index{typeid(Input)}))
            {
                ERR("Input {} of model {} has not been computed/set before the "
                    "model evaluation.",
                    BaseLib::typeToString<Input>(),
                    BaseLib::typeToString<Model>());
                all_inputs_computed = false;
            }
        });
    if (!all_inputs_computed)
    {
        return false;
    }

    using ModelOutputsWrapped = mp_transform<mp_identity, ModelOutputs>;

    // All outputs are "computed data", now.
    bool no_output_precomputed = true;
    mp_for_each<ModelOutputsWrapped>(
        [&computed_data,
         &no_output_precomputed]<typename Output>(mp_identity<Output>)
        {
            auto const [it, emplaced] = computed_data.emplace(typeid(Output));

            if (!emplaced)
            {
                ERR("Output {} of model {} is computed more than once.",
                    BaseLib::typeToString<Output>(),
                    BaseLib::typeToString<Model>());
                no_output_precomputed = false;
            }
        });

    return no_output_precomputed;
}

template <typename... Models>
bool isEvalOrderCorrectRT(boost::mp11::mp_list<Models...>,
                          std::unordered_set<std::type_index>&& computed_data)
{
    return (isEvalOrderCorrectRT<Models>(computed_data) && ...);
}
}  // namespace detail

/// Checks at runtime if the given \c Models are evaluated in the right order if
/// evaluated in the order in which they appear in the list of \c Models.
///
/// I.e., all input data of a model must have been computed before that model
/// will be evaluated and no two models must compute the same data.
///
/// The passed \c Inputs are data that already have been computed before the
/// first model is evaluated.
template <typename Models, typename Inputs>
bool isEvalOrderCorrectRT()  // RT for runtime
{
    using namespace boost::mp11;

    static_assert(mp_is_list<Models>::value);
    static_assert(mp_is_list<Inputs>::value);

    // Wrap inputs. The elements of InputsWrapped are default constructible.
    using InputsWrapped = mp_transform<mp_identity, Inputs>;

    // "Holds" all data that has been computed successively by the invoked
    // models.
    std::unordered_set<std::type_index> computed_data;

    // All inputs are considered "computed data".
    mp_for_each<InputsWrapped>(
        [&computed_data]<typename Input>(mp_identity<Input>)
        { computed_data.emplace(typeid(Input)); });

    return detail::isEvalOrderCorrectRT(mp_rename<Models, mp_list>{},
                                        std::move(computed_data));
}
}  // namespace ProcessLib::Graph
