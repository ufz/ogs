/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <boost/property_tree/ptree_fwd.hpp>

#pragma once

#define ASSERT_ARRAY_NEAR(E,A,N,eps)\
    for (std::size_t i=0; i<(unsigned)(N); i++) \
        ASSERT_NEAR((E)[i], (A)[i], (eps));

#define ASSERT_ARRAY_EQ(E,A,N)\
    for (std::size_t i=0; i<(unsigned)(N); i++) \
        ASSERT_EQ((E)[i], (A)[i]);

boost::property_tree::ptree
readXml(const char xml[]);
