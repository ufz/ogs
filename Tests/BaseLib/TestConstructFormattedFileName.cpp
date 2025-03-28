/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <string>

#include "BaseLib/FileTools.h"

TEST(BaseLib, constructFormattedFileName)
{
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "test_{:timestep}", "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("test_2", formatted_filename);
    }
    {
        auto const formatted_filename_time =
            BaseLib::constructFormattedFileName("test_{:0.5time}", "mesh_name",
                                                2, 0.2, 3, true);
        ASSERT_EQ("test_0.20000", formatted_filename_time);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "test_{:iteration}", "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("test_3", formatted_filename);
    }
    {
        auto const formatted_filename_timestep_time =
            BaseLib::constructFormattedFileName("test_{:timestep}_{:0.5time}",
                                                "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("test_2_0.20000", formatted_filename_timestep_time);
    }
    {
        auto const formatted_filename_time_timestep =
            BaseLib::constructFormattedFileName("test_{:.4time}_{:timestep}",
                                                "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("test_0.2000_2", formatted_filename_time_timestep);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "test_{:.4time}_{:timestep}", "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("test_0.2000_2", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "_ts_{:timestep}_t_{:.4time}", "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("_ts_2_t_0.2000", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "_ts_{:0>3timestep}_t_{:.4time}", "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("_ts_002_t_0.2000", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "_ts_{:0>3timestep}_t_{:.4etime}", "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("_ts_002_t_2.0000e-01", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "{:meshname}_ts_{:0>3timestep}_t_{:.4etime}_iter_{:0>2iteration}",
            "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("mesh_name_ts_002_t_2.0000e-01_iter_03", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "{:meshname}{:converged}", "mesh_name", 2, 0.2, 3, true);
        ASSERT_EQ("mesh_name", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "{:meshname}{:converged}", "mesh_name", 2, 0.2, 3, false);
        ASSERT_EQ("mesh_name_not_converged", formatted_filename);
    }
}
