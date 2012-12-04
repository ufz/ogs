/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file TestBidirectionalMap.cpp
 *
 * Created on 2012-11-03 by Norihiro Watanabe
 */

// ** INCLUDES **
#include "gtest/gtest.h"

#include <string>
#include "BidirectionalMap.h"

//-------------------------------------------------------------------
// Common testing functions for BidirectionalMap
//-------------------------------------------------------------------
// check mapping if the given key is not found in A
template <typename T1, typename T2>
static void checkNomatchingA(const BaseLib::BidirectionalMap<T1, T2> &map, const T1 &a)
{
    bool isExceptionThrowed = false;
    try {
        map.mapAtoB(a);
    } catch (std::exception &) {
        isExceptionThrowed = true;
    }
    ASSERT_TRUE(isExceptionThrowed);
}

// check mapping if the given key is not found in B
template <typename T1, typename T2>
static void checkNomatchingB(const BaseLib::BidirectionalMap<T1, T2> &map, const T2 &b)
{
    bool isExceptionThrowed = false;
    try {
        map.mapBtoA(b);
    } catch (std::out_of_range &) {
        isExceptionThrowed = true;
    }
    ASSERT_TRUE(isExceptionThrowed);
}

// simple testing assuming both keys are the same value but with different signs
static void checkIntIntState(const BaseLib::BidirectionalMap<int, int> &map, int i)
{
    const int j = -i;
    ASSERT_EQ(1, map.countInA(i));
    ASSERT_EQ(0, map.countInA(j));
    ASSERT_EQ(0, map.countInB(i));
    ASSERT_EQ(1, map.countInB(j));
    ASSERT_EQ(j, map.mapAtoB(i));
    checkNomatchingA(map, j);
    checkNomatchingB(map, i);
    ASSERT_EQ(i, map.mapBtoA(j));
};

static void checkStrStrState(const BaseLib::BidirectionalMap<std::string, std::string> &map, const std::string &str)
{
    const std::string str2 = "-" + str;
    ASSERT_EQ(1, map.countInA(str));
    ASSERT_EQ(0, map.countInA(str2));
    ASSERT_EQ(0, map.countInB(str));
    ASSERT_EQ(1, map.countInB(str2));
    ASSERT_STREQ(str2.c_str(), map.mapAtoB(str).c_str());
    checkNomatchingA(map, str2);
    checkNomatchingB(map, str);
    ASSERT_STREQ(str.c_str(), map.mapBtoA(str2).c_str());
};

static void checkIntStrState(const BaseLib::BidirectionalMap<int, std::string> &map, int i)
{
    std::stringstream ss;
    ss << -i;
    const std::string j = ss.str();

    ASSERT_EQ(1, map.countInA(i));
    ASSERT_EQ(0, map.countInA(-i));
    //ASSERT_EQ(0, map.countInB(i));
    ASSERT_EQ(1, map.countInB(j));
    ASSERT_EQ(j, map.mapAtoB(i));
    checkNomatchingA(map, -i);
    //checkNomatchingB(map, i);
    ASSERT_EQ(i, map.mapBtoA(j));
};

//-------------------------------------------------------------------
// Test cases
//-------------------------------------------------------------------
TEST(BaseLib, BidirectionalMap_IntInt_Normal) 
{
    BaseLib::BidirectionalMap<int, int> map1;

    // check initial
    ASSERT_EQ(0, map1.size());
    ASSERT_EQ(0, map1.countInA(1));
    ASSERT_EQ(0, map1.countInB(1));
    checkNomatchingA(map1, 1);
    checkNomatchingB(map1, 1);

    // should work without entries
    map1.clear();

    // check insert
    map1.insert(1, -1);
    ASSERT_EQ(1, map1.size());
    checkIntIntState(map1, 1);

    // check append
    map1.insert(2, -2);
    ASSERT_EQ(2, map1.size());
    checkIntIntState(map1, 1);
    checkIntIntState(map1, 2);

    // should work with entries
    map1.clear();
    ASSERT_EQ(0, map1.size());

    // check re-insert
    map1.insert(1, -1);
    ASSERT_EQ(1, map1.size());
    checkIntIntState(map1, 1);

    // check copy
    BaseLib::BidirectionalMap<int, int> map2(map1);
    ASSERT_EQ(1, map2.size());
    checkIntIntState(map2, 1);
    map2.clear();
    ASSERT_EQ(0, map2.size());
    ASSERT_EQ(1, map1.size());

    // check assignment
    BaseLib::BidirectionalMap<int, int> map3;
    map3 = map1;
    ASSERT_EQ(1, map3.size());
    checkIntIntState(map3, 1);
    map3.clear();
    ASSERT_EQ(0, map3.size());
    ASSERT_EQ(1, map1.size());

}

TEST(BaseLib, BidirectionalMap_StrStr_Normal) 
{
    BaseLib::BidirectionalMap<std::string, std::string> map1;

    // check initial
    ASSERT_EQ(0, map1.size());
    ASSERT_EQ(0, map1.countInA("1"));
    ASSERT_EQ(0, map1.countInB("1"));
    checkNomatchingA<std::string, std::string>(map1, "1");
    checkNomatchingB<std::string, std::string>(map1, "1");

    // should work without entries
    map1.clear();

    // check insert
    map1.insert("1", "-1");
    ASSERT_EQ(1, map1.size());
    checkStrStrState(map1, "1");

    // check append
    map1.insert("2", "-2");
    ASSERT_EQ(2, map1.size());
    checkStrStrState(map1, "1");
    checkStrStrState(map1, "2");

    // should work with entries
    map1.clear();
    ASSERT_EQ(0, map1.size());

    // check re-insert
    map1.insert("1", "-1");
    ASSERT_EQ(1, map1.size());
    checkStrStrState(map1, "1");

    // check copy
    BaseLib::BidirectionalMap<std::string, std::string> map2(map1);
    ASSERT_EQ(1, map2.size());
    checkStrStrState(map2, "1");
    map2.clear();
    ASSERT_EQ(0, map2.size());
    ASSERT_EQ(1, map1.size());

    // check assignment
    BaseLib::BidirectionalMap<std::string, std::string> map3;
    map3 = map1;
    ASSERT_EQ(1, map3.size());
    checkStrStrState(map3, "1");
    map3.clear();
    ASSERT_EQ(0, map3.size());
    ASSERT_EQ(1, map1.size());

    map1.clear();
}

TEST(BaseLib, BidirectionalMap_IntStr_Normal) 
{
    BaseLib::BidirectionalMap<int, std::string> map1;

    // check initial
    ASSERT_EQ(0, map1.size());
    ASSERT_EQ(0, map1.countInA(1));
    ASSERT_EQ(0, map1.countInB("1"));
    checkNomatchingA(map1, 1);
    checkNomatchingB<int, std::string>(map1, "1");

    // should work without entries
    map1.clear();

    // check insert
    map1.insert(1, "-1");
    ASSERT_EQ(1, map1.size());
    checkIntStrState(map1, 1);

    // check append
    map1.insert(2, "-2");
    ASSERT_EQ(2, map1.size());
    checkIntStrState(map1, 1);
    checkIntStrState(map1, 2);

    // should work with entries
    map1.clear();
    ASSERT_EQ(0, map1.size());

    // check re-insert
    map1.insert(1, "-1");
    ASSERT_EQ(1, map1.size());
    checkIntStrState(map1, 1);

    // check copy
    BaseLib::BidirectionalMap<int, std::string> map2(map1);
    ASSERT_EQ(1, map2.size());
    checkIntStrState(map2, 1);
    map2.clear();
    ASSERT_EQ(0, map2.size());
    ASSERT_EQ(1, map1.size());

    // check assignment
    BaseLib::BidirectionalMap<int, std::string> map3;
    map3 = map1;
    ASSERT_EQ(1, map3.size());
    checkIntStrState(map3, 1);
    map3.clear();
    ASSERT_EQ(0, map3.size());
    ASSERT_EQ(1, map1.size());

    map1.clear();
}
