/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <array>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <string>

#include <gtest/gtest.h>

#include "BaseLib/FileTools.h"

TEST(BaseLib, constructFormattedFileName)
{
    auto const formatted_filename =
        BaseLib::constructFormattedFileName("test_{:timestep}", 0, 2, 0.2);
    ASSERT_EQ("test_2", formatted_filename);

    auto const formatted_filename_time =
        BaseLib::constructFormattedFileName("test_{:0.5time}", 0, 2, 0.2);
    ASSERT_EQ("test_0.20000", formatted_filename_time);
    auto const formatted_filename_timestep_time =
        BaseLib::constructFormattedFileName("test_{:timestep}_{:0.5time}", 0, 2,
                                            0.2);
    ASSERT_EQ("test_2_0.20000", formatted_filename_timestep_time);
    auto const formatted_filename_time_timestep =
        BaseLib::constructFormattedFileName("test_{:.4time}_{:timestep}", 0, 2,
                                            0.2);
    ASSERT_EQ("test_0.2000_2", formatted_filename_time_timestep);

    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "test_{:.4time}_{:timestep}_{:process_id}", 0, 2, 0.2);
        ASSERT_EQ("test_0.2000_2_0", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "pcs_{:process_id}_ts_{:timestep}_t_{:.4time}", 0, 2, 0.2);
        ASSERT_EQ("pcs_0_ts_2_t_0.2000", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "pcs_{:process_id}_ts_{:0>3timestep}_t_{:.4time}", 0, 2, 0.2);
        ASSERT_EQ("pcs_0_ts_002_t_0.2000", formatted_filename);
    }
    {
        auto const formatted_filename = BaseLib::constructFormattedFileName(
            "pcs_{:process_id}_ts_{:0>3timestep}_t_{:.4etime}", 0, 2, 0.2);
        ASSERT_EQ("pcs_0_ts_002_t_2.0000e-01", formatted_filename);
    }
}
