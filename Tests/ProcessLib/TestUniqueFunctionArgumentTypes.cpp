/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Invoke.h"

TEST(ProcessLib_UniqueFunctionArgumentTypes, Unique)
{
    struct S
    {
        void no_args() {}              // support for zero argument methods
        void no_args_2() const {}      // support for const methods
        int no_args_3() { return 1; }  // support for arbitrary return types
        int& no_args_4() const
        {
            static int i;
            return i;
        }

        void one_arg(int) {}  // support for one argument methods
        double* one_arg_2(double& d) const
        {
            return &d;
        }  // argument type does not matter
        void one_arg_3(std::string const&) const {}  // classes work, too

        double two_args(int, double) { return 1.5; }
        int two_args_2(double*, int) const { return 1; }
        int two_args_3(int*, int) const
        {
            return 1;
        }  // difference in pointer/non-pointer is considered unique

        int three_args(std::string, int, char) const { return 1; }
        std::string four_args(unsigned, std::string&, double, char)
        {
            return "";
        }

        // To check the convenience function and macro
        void eval(unsigned, std::string const&, double*, char) {}
    };

    namespace D = ProcessLib::ThermoRichardsMechanics::detail;

    static_assert(D::areEvalArgumentTypesUnique(&S::no_args));
    static_assert(D::areEvalArgumentTypesUnique(&S::no_args_2));
    static_assert(D::areEvalArgumentTypesUnique(&S::no_args_3));
    static_assert(D::areEvalArgumentTypesUnique(&S::no_args_4));

    static_assert(D::areEvalArgumentTypesUnique(&S::one_arg));
    static_assert(D::areEvalArgumentTypesUnique(&S::one_arg_2));
    static_assert(D::areEvalArgumentTypesUnique(&S::one_arg_3));

    static_assert(D::areEvalArgumentTypesUnique(&S::two_args));
    static_assert(D::areEvalArgumentTypesUnique(&S::two_args_2));
    static_assert(D::areEvalArgumentTypesUnique(&S::two_args_3));

    static_assert(D::areEvalArgumentTypesUnique(&S::three_args));
    static_assert(D::areEvalArgumentTypesUnique(&S::four_args));

    static_assert(
        ProcessLib::ThermoRichardsMechanics::areEvalArgumentTypesUnique<S>());

    S const s;
    ProcessLib::ThermoRichardsMechanics::assertEvalArgsUnique(s);
}

TEST(ProcessLib_UniqueFunctionArgumentTypes, NonUnique)
{
    struct S
    {
        void same_type(int, int) {}

        void diff_in_const(double, const double) const {}
        void diff_in_ref(std::string, std::string&) const {}
        void diff_in_rvalue_ref(std::string, std::string&&) {}
        void diff_in_rvalue_ref_2(std::string const&, std::string&&) const {}
        void diff_in_const_ref(char, const char&) {}

        // same with other arguments interspersed
        void v2_same_type(char, int*, int*) {}

        void v2_diff_in_const(double*, int, double* const) const {}
        void v2_diff_in_ref(char, unsigned, std::string, std::string&) {}
        void v2_diff_in_rvalue_ref(std::string, std::string&&, double) const {}
        void v2_diff_in_rvalue_ref_2(long, std::string const&, std::string&&,
                                     double) const
        {
        }
        void v2_diff_in_const_ref(std::string&, char, int*, const char&,
                                  unsigned)
        {
        }

        // To check the convenience function
        void eval(unsigned, std::string const&, double*, char, std::string) {}
    };

    namespace D = ProcessLib::ThermoRichardsMechanics::detail;

    static_assert(!D::areEvalArgumentTypesUnique(&S::same_type));
    static_assert(!D::areEvalArgumentTypesUnique(&S::diff_in_const));
    static_assert(!D::areEvalArgumentTypesUnique(&S::diff_in_ref));
    static_assert(!D::areEvalArgumentTypesUnique(&S::diff_in_rvalue_ref));
    static_assert(!D::areEvalArgumentTypesUnique(&S::diff_in_rvalue_ref_2));
    static_assert(!D::areEvalArgumentTypesUnique(&S::diff_in_const_ref));

    static_assert(!D::areEvalArgumentTypesUnique(&S::v2_same_type));
    static_assert(!D::areEvalArgumentTypesUnique(&S::v2_diff_in_const));
    static_assert(!D::areEvalArgumentTypesUnique(&S::v2_diff_in_ref));
    static_assert(!D::areEvalArgumentTypesUnique(&S::v2_diff_in_rvalue_ref));
    static_assert(!D::areEvalArgumentTypesUnique(&S::v2_diff_in_rvalue_ref_2));
    static_assert(!D::areEvalArgumentTypesUnique(&S::v2_diff_in_const_ref));

    static_assert(
        !ProcessLib::ThermoRichardsMechanics::areEvalArgumentTypesUnique<S>());
}
