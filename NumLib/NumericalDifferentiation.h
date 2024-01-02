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

#include <Eigen/Core>
#include <tuple>
#include <utility>

#include "BaseLib/StrongType.h"

namespace NumLib
{
using RelativeEpsilon = BaseLib::StrongType<double, struct RelativeEpsilonTag>;
using MinimumPerturbation =
    BaseLib::StrongType<double, struct MinimumPerturbationTag>;

namespace detail
{
template <typename T>
struct IsScalar : std::true_type
{
};

template <int N>
struct IsScalar<Eigen::Matrix<double, N, 1, Eigen::ColMajor, N, 1>>
    : std::false_type
{
};

template <std::size_t IndexInTuple, typename Tuple>
double getScalarOrVectorComponent(Tuple const& tuple, Eigen::Index component)
{
    auto const& value = std::get<IndexInTuple>(tuple);

    if constexpr (IsScalar<std::remove_cvref_t<decltype(value)>>::value)
    {
        return value;
    }
    else
    {
        return value[component];
    }
}

/// Central differences derivative of a function with respect to a single scalar
/// variable.
///
/// See also NumLib::CentralDifferencesStrategy.
struct ComputeDerivativeWrtOneScalar_CD
{
    template <typename Function, typename TupleOfArgs,
              typename PerturbationStrategy, std::size_t PerturbedArgIdx,
              std::size_t... AllArgIdcs>
    auto operator()(Function const& f, TupleOfArgs const& args,
                    PerturbationStrategy const& pert_strat,
                    std::integral_constant<std::size_t, PerturbedArgIdx>,
                    Eigen::Index const perturbed_arg_component,
                    std::index_sequence<AllArgIdcs...>) const
    {
        auto const value_plus = f(pert_strat.perturbIf(
            std::bool_constant<PerturbedArgIdx == AllArgIdcs>{},
            std::get<AllArgIdcs>(args), 1.0, perturbed_arg_component)...);

        auto const value_minus = f(pert_strat.perturbIf(
            std::bool_constant<PerturbedArgIdx == AllArgIdcs>{},
            std::get<AllArgIdcs>(args), -1.0, perturbed_arg_component)...);

        auto const pert = pert_strat.getPerturbation(
            getScalarOrVectorComponent<PerturbedArgIdx>(
                args, perturbed_arg_component));

        // decltype enforces evaluation of Eigen expressions
        decltype(value_plus) deriv = (value_plus - value_minus) / (2 * pert);

        return deriv;
    }
};

/// Forward differences derivative of a function with respect to a single scalar
/// variable.
///
/// See also NumLib::ForwardDifferencesStrategy.
template <typename Value>
struct ComputeDerivativeWrtOneScalar_FD
{
    explicit ComputeDerivativeWrtOneScalar_FD(Value&& unperturbed_value)
        : unperturbed_value_{std::move(unperturbed_value)}
    {
    }

    template <typename Function, typename TupleOfArgs,
              typename PerturbationStrategy, std::size_t PerturbedArgIdx,
              std::size_t... AllArgIdcs>
    Value operator()(Function const& f, TupleOfArgs const& args,
                     PerturbationStrategy const& pert_strat,
                     std::integral_constant<std::size_t, PerturbedArgIdx>,
                     Eigen::Index const perturbed_arg_component,
                     std::index_sequence<AllArgIdcs...>) const
    {
        auto const value_plus = f(pert_strat.perturbIf(
            std::bool_constant<PerturbedArgIdx == AllArgIdcs>{},
            std::get<AllArgIdcs>(args), 1.0, perturbed_arg_component)...);

        auto const pert = pert_strat.getPerturbation(
            detail::getScalarOrVectorComponent<PerturbedArgIdx>(
                args, perturbed_arg_component));

        return (value_plus - unperturbed_value_) / pert;
    }

private:
    Value unperturbed_value_;
};

/// Perturbs values according to a relative epsilon. The perturbation should not
/// be less than a certain minimum.
struct DefaultPerturbationStrategy
{
    DefaultPerturbationStrategy(RelativeEpsilon const& rel_eps,
                                MinimumPerturbation const& min_pert)
        : rel_eps_{*rel_eps}, min_pert_{*min_pert}
    {
    }

    double getPerturbation(double const value) const
    {
        auto const pert = std::abs(value) * rel_eps_;

        if (std::abs(pert) >= std::abs(min_pert_))
        {
            return pert;
        }

        return min_pert_;
    }

    template <typename T>
    static T const& perturbIf(std::false_type, T const& value,
                              double const /*plus_or_minus*/,
                              Eigen::Index /*comp*/)
    {
        return value;
    }

    double perturbIf(std::true_type, double value, double const plus_or_minus,
                     Eigen::Index /*comp*/) const
    {
        return value + plus_or_minus * getPerturbation(value);
    }

    template <int N>
    Eigen::Vector<double, N> perturbIf(
        std::true_type,
        Eigen::Matrix<double, N, 1, Eigen::ColMajor, N, 1> const& vec,
        double const plus_or_minus,
        Eigen::Index comp) const
    {
        Eigen::Vector<double, N> vec_pert = vec;
        vec_pert[comp] += plus_or_minus * getPerturbation(vec[comp]);
        return vec_pert;
    }

private:
    double rel_eps_;
    double min_pert_;
};
}  // namespace detail

/// Strategy to compute numerical derivatives with central differences.
///
/// For use with NumLib::NumericalDerivative.
struct CentralDifferencesStrategy
{
    template <typename Function, typename... Args>
    static detail::ComputeDerivativeWrtOneScalar_CD createDByDScalar(
        Function const& /*f*/, Args const&... /*args*/)
    {
        return {};
    }
};

/// Strategy to compute numerical derivatives with forward differences.
///
/// For use with NumLib::NumericalDerivative.
struct ForwardDifferencesStrategy
{
    template <typename Function, typename... Args>
    static auto createDByDScalar(Function const& f, Args const&... args)
    {
        return detail::ComputeDerivativeWrtOneScalar_FD{f(args...)};
    }
};

// TODO better call it NumericalDifferentiationAlgorithm?
/// Computes numerical derivatives of almost arbitrary function (objects) with
/// respect to scalar or vectorial (Eigen vectors) variables.
///
/// \tparam DerivativeStrategy compute central or forward differences
template <typename DerivativeStrategy>
struct NumericalDerivative
{
    NumericalDerivative(RelativeEpsilon const& rel_eps,
                        MinimumPerturbation const& min_pert)
        : pert_strat_{rel_eps, min_pert}
    {
    }

    template <typename Function, typename... Args>
    auto operator()(Function const& f, Args const&... args) const
    {
        auto const d_by_dScalar =
            DerivativeStrategy::createDByDScalar(f, args...);

        // TODO also return value from the function, not only the derivatives?
        return differentiate(f,
                             std::forward_as_tuple(args...),
                             d_by_dScalar,
                             std::make_index_sequence<sizeof...(Args)>{});
    }

private:
    template <typename Function, typename TupleOfArgs, typename DByDScalar,
              std::size_t... AllArgIdcs>
    auto differentiate(Function const& f, TupleOfArgs const& args,
                       DByDScalar const& d_by_dScalar,
                       std::index_sequence<AllArgIdcs...> all_arg_idcs) const
    {
        return std::tuple{differentiateWrtScalarOrVectorialArgument(
            detail::IsScalar<std::remove_cvref_t<
                std::tuple_element_t<AllArgIdcs, TupleOfArgs>>>{},
            f, args, d_by_dScalar,
            std::integral_constant<std::size_t, AllArgIdcs>{},
            all_arg_idcs)... /* "for each function argument" */};
    }

    // scalar case
    template <typename Function, typename TupleOfArgs, typename DByDScalar,
              std::size_t... AllArgIdcs, std::size_t PerturbedArgIdx>
    auto differentiateWrtScalarOrVectorialArgument(
        std::true_type /* is_scalar */, Function const& f,
        TupleOfArgs const& args, DByDScalar const& d_by_dScalar,
        std::integral_constant<std::size_t, PerturbedArgIdx> perturbed_arg_idx,
        std::index_sequence<AllArgIdcs...> all_arg_idcs) const
    {
        constexpr Eigen::Index component_does_not_matter = -1;

        return d_by_dScalar(f, args, pert_strat_, perturbed_arg_idx,
                            component_does_not_matter, all_arg_idcs);
    }

    // vectorial case
    template <typename Function, typename TupleOfArgs, typename DByDScalar,
              std::size_t... AllArgIdcs, std::size_t PerturbedArgIdx>
    auto differentiateWrtScalarOrVectorialArgument(
        std::false_type /* is_scalar */, Function const& f,
        TupleOfArgs const& args, DByDScalar const& d_by_dScalar,
        std::integral_constant<std::size_t, PerturbedArgIdx> perturbed_arg_idx,
        std::index_sequence<AllArgIdcs...> all_arg_idcs) const
    {
        using VectorialArg = std::remove_cvref_t<
            std::tuple_element_t<PerturbedArgIdx, TupleOfArgs>>;
        constexpr int N = VectorialArg::RowsAtCompileTime;

        static_assert(N != Eigen::Dynamic);
        static_assert(VectorialArg::ColsAtCompileTime == 1,
                      "Row vectors are not supported, yet. If you implement "
                      "support for them, make sure to test your implementation "
                      "thoroughly.");

        return differentiateWrtAllVectorComponents(
            f, args, d_by_dScalar,
            std::make_integer_sequence<Eigen::Index, N>{}, perturbed_arg_idx,
            all_arg_idcs);
    }

    template <typename Function, typename TupleOfArgs, typename DByDScalar,
              Eigen::Index... PerturbedArgComponents, std::size_t... AllArgIdcs,
              std::size_t PerturbedArgIdx>
    auto differentiateWrtAllVectorComponents(
        Function const& f, TupleOfArgs const& args,
        DByDScalar const& d_by_dScalar,
        std::integer_sequence<Eigen::Index, PerturbedArgComponents...>,
        std::integral_constant<std::size_t, PerturbedArgIdx> perturbed_arg_idx,
        std::index_sequence<AllArgIdcs...> all_arg_idcs) const
    {
        return std::array{
            d_by_dScalar(f, args, pert_strat_, perturbed_arg_idx,
                         PerturbedArgComponents, all_arg_idcs)...
            /* "for each component of the vectorial function argument being
               perturbed" */
        };
    }

    detail::DefaultPerturbationStrategy pert_strat_;
};

}  // namespace NumLib
