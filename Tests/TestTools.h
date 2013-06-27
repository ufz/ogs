/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#ifndef TESTTOOLS_H_
#define TESTTOOLS_H_

inline void ASSERT_DOUBLE_ARRAY_EQ(const double* Expected, const double* Actual, std::size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon);
}

template <typename T1, typename T2>
inline void ASSERT_DOUBLE_ARRAY_EQ(const T1 &Expected, const T2 &Actual, std::size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon);
}

template <typename T>
class series
{
public:
    series(T v0=0, T by=1) : v(v0), dv(by) {};
    T operator()() { v+=dv; return (v-dv);}
private:
    T v;
    T dv;
};
#endif // TESTTOOLS_H_
