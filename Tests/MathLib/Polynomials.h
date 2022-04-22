/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <gtest/gtest.h>

#include "Tests/VectorUtils.h"

namespace TestPolynomials
{

struct FBase
{
private:
    using Function = std::function<double(std::array<double, 3> const&)>;

    static std::vector<double> initCoeffs(std::size_t const num_coeffs)
    {
        std::vector<double> cs(num_coeffs);
        fillVectorRandomly(cs);
        return cs;
    }

protected:
    explicit FBase(std::size_t const num_coeffs)
        : coeffs(initCoeffs(num_coeffs))
    {
    }

public:
    virtual double operator()(
        std::array<double, 3> const& /*coords*/) const = 0;

    virtual double getAnalyticalIntegralOverUnitCube() const = 0;

    Function toFunction() const
    {
        // Note: this is captured, hence the caller is responsible that 'this'
        // is alive as long as the returned function is used.
        return [this](std::array<double, 3> const& coords)
        { return (*this)(coords); };
    }

    virtual ~FBase() = default;

    std::vector<double> const coeffs;
};

std::ostream& operator<<(std::ostream& os, FBase const& f)
{
    os << "Function{coeffs[" << f.coeffs.size() << "]={";
    for (std::size_t i = 0; i < f.coeffs.size(); ++i)
    {
        if (i != 0)
            os << ", ";
        os << f.coeffs[i];
    }
    os << "}}";
    return os;
}

// A constant function.
struct FConst final : FBase
{
    FConst() : FBase(1) {}

    double operator()(std::array<double, 3> const& /*unused*/) const override
    {
        return coeffs[0];
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        return coeffs[0];
    }
};

namespace detail
{
// Non-templated core implementation.
struct FLin : FBase
{
    FLin(std::size_t dim_) : FBase(1 + dim_), dim{dim_} {}

    double operator()(std::array<double, 3> const& coords) const override final
    {
        double linear_terms = 0.;

        for (std::size_t c = 0; c < dim && c < 3; ++c)
        {
            linear_terms += coeffs[c + 1] * coords[c];
        }
        return coeffs[0] + linear_terms;
    }

    double getAnalyticalIntegralOverUnitCube() const override final
    {
        return coeffs[0];
    }

private:
    std::size_t const dim;
};
}  // namespace detail

// A linear affine function f(x, y, ...) = c0 + c1 * x + c2 * y + ...
//
// Template redirects to non-templated core implementation.
template <std::size_t DIM>
#if __cpp_concepts >= 201907L
requires requires
{
    DIM <= 3;
}
#endif
struct FLin final : detail::FLin
{
    FLin() : detail::FLin{DIM} {}
};

// Quadratic function.
template <std::size_t DIM>
struct FQuad;

template <>
struct FQuad<3> final : FBase
{
    FQuad() : FBase(10) {}

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];
        auto const y = coords[1];
        auto const z = coords[2];
        return coeffs[0] + coeffs[1] * x + coeffs[2] * y + coeffs[3] * z +
               coeffs[4] * x * y + coeffs[5] * y * z + coeffs[6] * z * x +
               coeffs[7] * x * x + coeffs[8] * y * y + coeffs[9] * z * z;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double const a3 = a * a * a;
        double const b3 = b * b * b;

        return coeffs[0] + (coeffs[7] + coeffs[8] + coeffs[9]) * (b3 - a3) / 3.;
    }
};

template <>
struct FQuad<2> final : FBase
{
    FQuad() : FBase(6) {}

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];
        auto const y = coords[1];
        return coeffs[0] + coeffs[1] * x + coeffs[2] * y + coeffs[3] * x * y +
               coeffs[4] * x * x + coeffs[5] * y * y;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double const a3 = a * a * a;
        double const b3 = b * b * b;

        return coeffs[0] + (coeffs[4] + coeffs[5]) * (b3 - a3) / 3.;
    }
};

template <>
struct FQuad<1> final : FBase
{
    FQuad() : FBase(3) {}

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];
        return coeffs[0] + coeffs[1] * x + coeffs[2] * x * x;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double const a3 = a * a * a;
        double const b3 = b * b * b;

        return coeffs[0] + coeffs[2] * (b3 - a3) / 3.;
    }
};

namespace detail
{
// Non-templated core implementation.
struct FSeparablePolynomial : FBase
{
    explicit FSeparablePolynomial(unsigned dim, unsigned polynomial_degree)
        : FBase(dim * polynomial_degree + dim),
          _dim{dim},
          _degree(polynomial_degree)
    {
    }

    // f(x, y, z) = g(x) * h(y) * i(z)
    double operator()(std::array<double, 3> const& coords) const override final
    {
        double res = 1.0;
        for (unsigned d = 0; d < _dim; ++d)
        {
            auto const x = coords[d];

            double poly = 0.0;
            for (unsigned n = 0; n <= _degree; ++n)
            {
                poly += coeffs[n + d * (_degree + 1)] * std::pow(x, n);
            }

            res *= poly;
        }

        return res;
    }

    // [ F(x, y, z) ]_a^b = [ G(x) ]_a^b * [ H(y) ]_a^b * [ I(z) ]_a^b
    double getAnalyticalIntegralOverUnitCube() const override final
    {
        double const a = -.5;
        double const b = .5;

        double res = 1.0;
        for (unsigned d = 0; d < _dim; ++d)
        {
            double poly = 0.0;
            for (unsigned n = 0; n <= _degree; ++n)
            {
                poly += coeffs[n + d * (_degree + 1)] *
                        (std::pow(b, n + 1) - std::pow(a, n + 1)) / (n + 1);
            }

            res *= poly;
        }

        return res;
    }

private:
    unsigned const _dim;
    unsigned const _degree;
};

}  // namespace detail

// A polynomial in DIM variables that can be decomposed into a product of DIM
// polynomials in 1 variable: f(x, y, z) = g(x) * h(y) * i(z)
//
// Template redirects to non-templated core implementation.
template <std::size_t DIM>
#if __cpp_concepts >= 201907L
requires requires
{
    DIM <= 3;
}
#endif
struct FSeparablePolynomial final : detail::FSeparablePolynomial
{
    explicit FSeparablePolynomial(unsigned polynomial_degree)
        : detail::FSeparablePolynomial(DIM, polynomial_degree)
    {
    }
};

unsigned binomial_coefficient(unsigned n, unsigned k)
{
    EXPECT_GE(n, k);
    unsigned res = 1;

    for (unsigned i = n; i > k; --i)
    {
        res *= i;
    }

    for (unsigned i = n - k; i > 0; --i)
    {
        res /= i;
    }

    return res;
}

/* This function is a polynomial where for each monomial a_ijk x^i y^j z^k
 * holds: i + j + k <= n, where n is the overall polynomial degree
 */
template <std::size_t DIM>
struct FNonseparablePolynomial;

template <>
struct FNonseparablePolynomial<3> final : FBase
{
    // The number of coefficients/monomials are obtained as follows: Compute the
    // number of combinations with repititions when drawing
    // polynomial_degree times from the set { x, y, z, 1 }
    explicit FNonseparablePolynomial(unsigned polynomial_degree)
        : FBase(binomial_coefficient(4 + polynomial_degree - 1, 4 - 1)),
          _degree(polynomial_degree)
    {
    }

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];
        auto const y = coords[1];
        auto const z = coords[2];

        double res = 0.0;
        unsigned index = 0;

        for (unsigned x_deg = 0; x_deg <= _degree; ++x_deg)
        {
            for (unsigned y_deg = 0; x_deg + y_deg <= _degree; ++y_deg)
            {
                for (unsigned z_deg = 0; x_deg + y_deg + z_deg <= _degree;
                     ++z_deg)
                {
                    EXPECT_GT(coeffs.size(), index);

                    res += coeffs[index] * std::pow(x, x_deg) *
                           std::pow(y, y_deg) * std::pow(z, z_deg);

                    ++index;
                }
            }
        }

        EXPECT_EQ(coeffs.size(), index);

        return res;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double res = 0.0;
        unsigned index = 0;

        for (unsigned x_deg = 0; x_deg <= _degree; ++x_deg)
        {
            for (unsigned y_deg = 0; x_deg + y_deg <= _degree; ++y_deg)
            {
                for (unsigned z_deg = 0; x_deg + y_deg + z_deg <= _degree;
                     ++z_deg)
                {
                    EXPECT_GT(coeffs.size(), index);

                    res += coeffs[index] *
                           (std::pow(b, x_deg + 1) - std::pow(a, x_deg + 1)) /
                           (x_deg + 1) *
                           (std::pow(b, y_deg + 1) - std::pow(a, y_deg + 1)) /
                           (y_deg + 1) *
                           (std::pow(b, z_deg + 1) - std::pow(a, z_deg + 1)) /
                           (z_deg + 1);

                    ++index;
                }
            }
        }

        EXPECT_EQ(coeffs.size(), index);

        return res;
    }

private:
    unsigned const _degree;
};

template <>
struct FNonseparablePolynomial<2> final : FBase
{
    // The number of coefficients/monomials are obtained as follows: Compute the
    // number of combinations with repititions when drawing
    // polynomial_degree times from the set { x, y, 1 }
    explicit FNonseparablePolynomial(unsigned polynomial_degree)
        : FBase(binomial_coefficient(3 + polynomial_degree - 1, 3 - 1)),
          _degree(polynomial_degree)
    {
    }

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];
        auto const y = coords[1];

        double res = 0.0;
        unsigned index = 0;

        for (unsigned x_deg = 0; x_deg <= _degree; ++x_deg)
        {
            for (unsigned y_deg = 0; x_deg + y_deg <= _degree; ++y_deg)
            {
                EXPECT_GT(coeffs.size(), index);

                res += coeffs[index] * std::pow(x, x_deg) * std::pow(y, y_deg);

                ++index;
            }
        }

        EXPECT_EQ(coeffs.size(), index);

        return res;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double res = 0.0;
        unsigned index = 0;

        for (unsigned x_deg = 0; x_deg <= _degree; ++x_deg)
        {
            for (unsigned y_deg = 0; x_deg + y_deg <= _degree; ++y_deg)
            {
                EXPECT_GT(coeffs.size(), index);

                res += coeffs[index] *
                       (std::pow(b, x_deg + 1) - std::pow(a, x_deg + 1)) /
                       (x_deg + 1) *
                       (std::pow(b, y_deg + 1) - std::pow(a, y_deg + 1)) /
                       (y_deg + 1);

                ++index;
            }
        }

        EXPECT_EQ(coeffs.size(), index);

        return res;
    }

private:
    unsigned const _degree;
};

template <>
struct FNonseparablePolynomial<1> final : FBase
{
    // The number of coefficients/monomials are obtained as follows: Compute the
    // number of combinations with repititions when drawing
    // polynomial_degree times from the set { x, 1 }
    explicit FNonseparablePolynomial(unsigned polynomial_degree)
        : FBase(binomial_coefficient(2 + polynomial_degree - 1, 2 - 1)),
          _degree(polynomial_degree)
    {
    }

    double operator()(std::array<double, 3> const& coords) const override
    {
        auto const x = coords[0];

        double res = 0.0;
        unsigned index = 0;

        for (unsigned x_deg = 0; x_deg <= _degree; ++x_deg)
        {
            EXPECT_GT(coeffs.size(), index);

            res += coeffs[index] * std::pow(x, x_deg);

            ++index;
        }

        EXPECT_EQ(coeffs.size(), index);

        return res;
    }

    double getAnalyticalIntegralOverUnitCube() const override
    {
        double const a = -.5;
        double const b = .5;

        double res = 0.0;
        unsigned index = 0;

        for (unsigned x_deg = 0; x_deg <= _degree; ++x_deg)
        {
            EXPECT_GT(coeffs.size(), index);

            res += coeffs[index] *
                   (std::pow(b, x_deg + 1) - std::pow(a, x_deg + 1)) /
                   (x_deg + 1);

            ++index;
        }

        EXPECT_EQ(coeffs.size(), index);

        return res;
    }

private:
    unsigned const _degree;
};

template <std::size_t dim>
std::unique_ptr<FBase> getF(unsigned polynomial_order)
{
    switch (polynomial_order)
    {
        case 0:
            return std::make_unique<FConst>();
        case 1:
            return std::make_unique<FLin<dim>>();
        case 2:
            return std::make_unique<FQuad<dim>>();
    }

    OGS_FATAL("unsupported polynomial order: {:d}.", polynomial_order);
}

}  // namespace TestPolynomials
