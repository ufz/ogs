/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "BaseLib/StrongType.h"
#include "ProcessLib/Graph/Get.h"

TEST(ProcessLib_Graph_Get, ReadSingleTuple)
{
    namespace PG = ProcessLib::Graph;
    std::tuple<char, double, int> const t{'a', 5.5, 3};

    // check same values
    EXPECT_EQ('a', PG::get<char>(t));
    EXPECT_EQ(5.5, PG::get<double>(t));
    EXPECT_EQ(3, PG::get<int>(t));

    // check same addresses
    EXPECT_EQ(&std::get<char>(t), &PG::get<char>(t));
    EXPECT_EQ(&std::get<double>(t), &PG::get<double>(t));
    EXPECT_EQ(&std::get<int>(t), &PG::get<int>(t));
}

TEST(ProcessLib_Graph_Get, WriteSingleTuple)
{
    namespace PG = ProcessLib::Graph;
    std::tuple<char, double, int> t{'a', 5.5, 3};

    PG::get<char>(t) = 'b';
    PG::get<double>(t) = 6.25;
    PG::get<int>(t) = 5;

    // check with std::get
    EXPECT_EQ('b', std::get<char>(t));
    EXPECT_EQ(6.25, std::get<double>(t));
    EXPECT_EQ(5, std::get<int>(t));
}

TEST(ProcessLib_Graph_Get, ReadTwoTuples)
{
    namespace PG = ProcessLib::Graph;

    using S1 = BaseLib::StrongType<int, struct S1Tag>;
    using S2 = BaseLib::StrongType<double, struct S2Tag>;

    std::tuple<char, S1> const t1{'a', 2};
    std::tuple<S2, int> const t2{5.5, 4};

    // check same values
    EXPECT_EQ('a', PG::get<char>(t1, t2));
    EXPECT_EQ(2, *PG::get<S1>(t1, t2));
    EXPECT_EQ(5.5, *PG::get<S2>(t1, t2));
    EXPECT_EQ(4, PG::get<int>(t1, t2));

    // order does not matter
    EXPECT_EQ('a', PG::get<char>(t2, t1));
    EXPECT_EQ(2, *PG::get<S1>(t2, t1));
    EXPECT_EQ(5.5, *PG::get<S2>(t2, t1));
    EXPECT_EQ(4, PG::get<int>(t2, t1));

    // check same addresses
    EXPECT_EQ(&std::get<char>(t1), &PG::get<char>(t1, t2));
    EXPECT_EQ(&std::get<S1>(t1), &PG::get<S1>(t1, t2));
    EXPECT_EQ(&std::get<S2>(t2), &PG::get<S2>(t1, t2));
    EXPECT_EQ(&std::get<int>(t2), &PG::get<int>(t1, t2));
}

TEST(ProcessLib_Graph_Get, WriteTwoTuples)
{
    namespace PG = ProcessLib::Graph;

    using S1 = BaseLib::StrongType<int, struct S1Tag>;
    using S2 = BaseLib::StrongType<double, struct S2Tag>;

    std::tuple<char, S1> t1{'a', 2};
    std::tuple<S2, int> t2{5.5, 4};

    PG::get<char>(t1, t2) = 'c';
    PG::get<S1>(t2, t1) = S1{7};
    PG::get<S2>(t1, t2) = S2{6.5};
    PG::get<int>(t2, t1) = 9;

    // check with std::get
    EXPECT_EQ('c', std::get<char>(t1));
    EXPECT_EQ(7, *std::get<S1>(t1));
    EXPECT_EQ(6.5, *std::get<S2>(t2));
    EXPECT_EQ(9, std::get<int>(t2));
}

TEST(ProcessLib_Graph_Get, ReadMultipleTuples)
{
    namespace PG = ProcessLib::Graph;

    using S1 = BaseLib::StrongType<int, struct S1Tag>;
    using S2 = BaseLib::StrongType<double, struct S2Tag>;
    using S3 = BaseLib::StrongType<char, struct S3Tag>;

    std::tuple<char, S1> const t1{'a', 2};
    std::tuple<S2> const t2{5.5};
    std::tuple<> const t3;
    std::tuple<std::string, S3, int> const t4{"str", 'x', 4};

    // check same values
    EXPECT_EQ('a', PG::get<char>(t1, t2, t3, t4));
    EXPECT_EQ(2, *PG::get<S1>(t3, t1, t2, t4));
    EXPECT_EQ(5.5, *PG::get<S2>(t4, t1, t3, t2));
    EXPECT_EQ(4, PG::get<int>(t1, t4, t3, t2));
    EXPECT_EQ("str", PG::get<std::string>(t3, t4, t1, t2));
    EXPECT_EQ('x', *PG::get<S3>(t4, t3, t1, t2));

    // check same addresses
    EXPECT_EQ(&std::get<char>(t1), &PG::get<char>(t1, t2, t3, t4));
    EXPECT_EQ(&std::get<S1>(t1), &PG::get<S1>(t3, t4, t1, t2));
    EXPECT_EQ(&std::get<S2>(t2), &PG::get<S2>(t1, t3, t2, t4));
    EXPECT_EQ(&std::get<int>(t4), &PG::get<int>(t1, t4, t2, t3));
    EXPECT_EQ(&std::get<S3>(t4), &PG::get<S3>(t2, t4, t1, t3));
    EXPECT_EQ(&std::get<std::string>(t4),
              &PG::get<std::string>(t2, t3, t1, t4));
}

TEST(ProcessLib_Graph_Get, WriteMultipleTuples)
{
    namespace PG = ProcessLib::Graph;

    using S1 = BaseLib::StrongType<int, struct S1Tag>;
    using S2 = BaseLib::StrongType<double, struct S2Tag>;
    using S3 = BaseLib::StrongType<char, struct S3Tag>;

    std::tuple<char, S1> t1{'a', 2};
    std::tuple<S2> t2{5.5};
    std::tuple<> t3;
    std::tuple<std::string, S3, int> t4{"str", 'x', 4};

    PG::get<char>(t1, t2, t3, t4) = 'c';
    PG::get<S1>(t3, t1, t2, t4) = S1{5};
    PG::get<S2>(t4, t1, t3, t2) = S2{7.5};
    PG::get<int>(t1, t4, t3, t2) = 11;
    PG::get<std::string>(t3, t4, t1, t2) = "rts";
    PG::get<S3>(t4, t3, t1, t2) = S3{'y'};

    EXPECT_EQ('c', std::get<char>(t1));
    EXPECT_EQ(5, *std::get<S1>(t1));
    EXPECT_EQ(7.5, *std::get<S2>(t2));
    EXPECT_EQ(11, std::get<int>(t4));
    EXPECT_EQ("rts", std::get<std::string>(t4));
    EXPECT_EQ('y', *std::get<S3>(t4));
}

TEST(ProcessLib_Graph_Get, ReadTwoTuplesWithReferenceMembers)
{
    namespace PG = ProcessLib::Graph;

    using S1 = BaseLib::StrongType<int, struct S1Tag>;
    using S2 = BaseLib::StrongType<double, struct S2Tag>;

    S1 i{2};
    S2 const d{5.5};
    std::tuple<char, S1&> const t1{'a', i};
    std::tuple<S2 const&, int> const t2{d, 4};

    // check same values
    EXPECT_EQ('a', PG::get<char>(t1, t2));
    EXPECT_EQ(2, *PG::get<S1>(t2, t1));
    EXPECT_EQ(5.5, *PG::get<S2>(t1, t2));
    EXPECT_EQ(4, PG::get<int>(t2, t1));

    // check same addresses
    EXPECT_EQ(&std::get<char>(t1), &PG::get<char>(t1, t2));
    // note: PG::get<>() access without cvref in contrast to std::get<>().
    EXPECT_EQ(&std::get<S1&>(t1), &PG::get<S1>(t1, t2));
    EXPECT_EQ(&std::get<S2 const&>(t2), &PG::get<S2>(t1, t2));
    EXPECT_EQ(&std::get<int>(t2), &PG::get<int>(t1, t2));

    EXPECT_EQ(&i, &PG::get<S1>(t1, t2));
    EXPECT_EQ(&d, &PG::get<S2>(t1, t2));
}

TEST(ProcessLib_Graph_Get, WriteTwoTuplesWithReferenceMembers)
{
    namespace PG = ProcessLib::Graph;

    using S1 = BaseLib::StrongType<int, struct S1Tag>;
    using S2 = BaseLib::StrongType<double, struct S2Tag>;

    S1 i{2};
    S2 d{5.5};
    std::tuple<char, S1&> t1{'a', i};
    std::tuple<S2&, int> t2{d, 4};

    PG::get<char>(t1, t2) = 'c';
    PG::get<S1>(t2, t1) = S1{7};
    PG::get<S2>(t1, t2) = S2{6.5};
    PG::get<int>(t2, t1) = 9;

    EXPECT_EQ('c', std::get<char>(t1));
    EXPECT_EQ(7, *i);  // original modified
    EXPECT_EQ(6.5, *d);
    EXPECT_EQ(9, std::get<int>(t2));
}
