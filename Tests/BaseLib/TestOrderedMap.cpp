/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file TestOrderedMap.cpp
 *
 * Created on 2012-11-03 by Norihiro Watanabe
 */

// ** INCLUDES **
#include "gtest.h"

#include "OrderedMap.h"

//-------------------------------------------------------------------
// Common testing functions for BidirectionalMap
//-------------------------------------------------------------------
template <typename T1, typename T2>
static void checkOutOfRange(const BaseLib::OrderedMap<T1, T2> &map, const size_t i)
{
    bool isExceptionThrowed = false;
    try {
        map[i]->second;
    } catch (std::exception &) {
        isExceptionThrowed = true;
    }
    ASSERT_TRUE(isExceptionThrowed);
}

// simple testing assuming both keys are the same value but with different signs
static void checkIntIntState(const BaseLib::OrderedMap<int, int> &map, int i)
{
    const int j = -i;
    ASSERT_EQ(1, map.count(i));
    ASSERT_EQ(0, map.count(-i));
    ASSERT_EQ(j, map.find(i)->second);
};

//-------------------------------------------------------------------
// Test cases
//-------------------------------------------------------------------
TEST(BaseLib, OrderedMap_IntInt_Normal) 
{
    BaseLib::OrderedMap<int, int> map1;
    
    // check initial
    ASSERT_EQ(0, map1.size());
    ASSERT_EQ(0, map1.count(1));
    ASSERT_EQ(map1.end(), map1.find(1));
    checkOutOfRange(map1, 1);

    // should work without entries
    map1.clear();

    // check insert
    map1.insert(2, -2);
    ASSERT_EQ(1, map1.size());
    checkIntIntState(map1, 2);
    ASSERT_EQ(2, map1[0]->first);

    // check append
    map1.insert(3, -3);
    map1.insert(1, -1);
    ASSERT_EQ(3, map1.size());
    checkIntIntState(map1, 2);
    checkIntIntState(map1, 3);
    checkIntIntState(map1, 1);
    ASSERT_EQ(2, map1[0]->first);
    ASSERT_EQ(3, map1[1]->first);
    ASSERT_EQ(1, map1[2]->first);

    // should work with entries
    map1.clear();
    ASSERT_EQ(0, map1.size());

    // check re-insert
    map1.insert(1, -1);
    ASSERT_EQ(1, map1.size());
    checkIntIntState(map1, 1);
    ASSERT_EQ(1, map1[0]->first);

}

