/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-12-05
 * \brief  Implementation of the LisOption class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisOption.h"

#include <logog/include/logog.hpp>

namespace MathLib
{

void LisOption::addOptions(const BaseLib::ConfigTree & options)
{
    auto const range = options.equal_range("lis_option");
    for (auto it=range.first; it!=range.second; ++it)
    {
        auto key = it->second.get_optional<Key>("option");
        auto val = it->second.get_optional<Value>("value");

        if (key && val) {
            settings[*key] = *val;
        }
    }

    auto solver_type = options.get_optional<std::string>("solver_type");
    if (solver_type) {
        settings["-i"] = *solver_type;
    }
    auto precon_type = options.get_optional<std::string>("precon_type");
    if (precon_type) {
        settings["-p"] = *precon_type;
    }
    auto error_tolerance = options.get_optional<std::string>("error_tolerance");
    if (error_tolerance) {
        settings["-tol"] = *error_tolerance;
    }
    auto max_iteration_step = options.get_optional<std::string>("max_iteration_step");
    if (max_iteration_step) {
        settings["-maxiter"] = *max_iteration_step;
    }
}

void LisOption::printInfo() const
{
    if (settings.empty()) {
        INFO("Lis options: none set.");
    } else {
        INFO("Lis options:");

        for (auto const& it : settings) {
            INFO("-> %s %s ", it.first.c_str(), it.second.c_str());
        }
        INFO("");
    }
}

} //MathLib
