/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <boost/property_tree/ptree_fwd.hpp>

#pragma once

#define ASSERT_ARRAY_NEAR(E, A, N, eps)             \
    for (std::size_t i = 0; i < (unsigned)(N); i++) \
        ASSERT_NEAR((E)[i], (A)[i], (eps));

#define ASSERT_ARRAY_EQ(E, A, N)                    \
    for (std::size_t i = 0; i < (unsigned)(N); i++) \
        ASSERT_EQ((E)[i], (A)[i]);

namespace Tests
{
boost::property_tree::ptree readXml(const char xml[]);

// A googletest predicate meant to be used with Eigen types with
// ASSERT_PRED_FORMAT2 and ASSERT_PRED_FORMAT3
struct EigenIsNear
{
    template <typename DerivedA, typename DerivedB>
    testing::AssertionResult operator()(const char* a_expr, const char* b_expr,
                                        const char* /*abstol_expr*/,
                                        Eigen::MatrixBase<DerivedA> const& a,
                                        Eigen::MatrixBase<DerivedB> const& b,
                                        double const abstol) const
    {
        static_assert(static_cast<Eigen::Index>(DerivedA::RowsAtCompileTime) ==
                          DerivedB::RowsAtCompileTime ||
                      DerivedA::RowsAtCompileTime == Eigen::Dynamic ||
                      DerivedB::RowsAtCompileTime == Eigen::Dynamic);
        static_assert(static_cast<Eigen::Index>(DerivedA::ColsAtCompileTime) ==
                          DerivedB::ColsAtCompileTime ||
                      DerivedA::ColsAtCompileTime == Eigen::Dynamic ||
                      DerivedB::ColsAtCompileTime == Eigen::Dynamic);

        if (a.rows() != b.rows())
        {
            return testing::AssertionFailure()
                   << a_expr << " and " << b_expr
                   << " have a different number of rows:  " << a.rows()
                   << " != " << b.rows() << '\n'
                   << a_expr << " evaluates to " << a << '\n'
                   << b_expr << " evaluates to " << b;
        }

        if (a.cols() != b.cols())
        {
            return testing::AssertionFailure()
                   << a_expr << " and " << b_expr
                   << " have a different number of columns:  " << a.cols()
                   << " != " << b.cols() << '\n'
                   << a_expr << " evaluates to " << a << '\n'
                   << b_expr << " evaluates to " << b;
        }

        auto const diff = (a - b).eval();

        for (Eigen::Index r = 0; r < a.rows(); ++r)
        {
            for (Eigen::Index c = 0; c < a.cols(); ++c)
            {
                auto const diff_comp = diff(r, c);

                if (!(std::abs(diff_comp) <= abstol)
                    // writing the comparison in this way also works with NaN
                )
                {
                    return testing::AssertionFailure()
                           << a_expr << " and " << b_expr << " differ by |"
                           << diff_comp << "| > " << abstol << " at position ("
                           << r << ", " << c << ")\n"
                           << a_expr << " evaluates to " << a << '\n'
                           << b_expr << " evaluates to " << b;
                }
            }
        }

        return testing::AssertionSuccess();
    }

    template <typename DerivedA, typename DerivedB>
    testing::AssertionResult operator()(
        const char* a_expr, const char* b_expr,
        Eigen::MatrixBase<DerivedA> const& a,
        Eigen::MatrixBase<DerivedB> const& b) const
    {
        return (*this)(a_expr, b_expr, "NO EXPR AVAIL", a, b,
                       std::numeric_limits<double>::epsilon());
    }
};
}  // namespace Tests
