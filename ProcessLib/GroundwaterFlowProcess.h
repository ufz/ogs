/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "MeshLib/Mesh.h"
#include "ProcessVariable.h"

namespace ProcessLib
{

class GroundwaterFlowProcess : public Process
{
    using ConfigTree = boost::property_tree::ptree;

public:
    GroundwaterFlowProcess(MeshLib::Mesh const& mesh,
            std::vector<ProcessVariable> const& variables,
            ConfigTree const& config)
        : Process(mesh, 1)
    {
        DBUG("Create GroundwaterFlowProcess.");

        // Find the corresponding process variable.
        std::string const name = config.get<std::string>("process_variable");

        auto const& variable = std::find_if(variables.cbegin(), variables.cend(),
                [&name](ProcessVariable const& v) {
                    return v.getName() == name;
                });

        if (variable == variables.end())
            ERR("Expected process variable \'%s\' not found in provided variables list.",
                name.c_str());

        DBUG("Associate hydraulic_head with process variable \'%s\'.",
            name.c_str());
        _hydraulic_head = &*variable;

    }

    void initialize()
    {
    }

    void solve()
    {
    }

    void post()
    {
    }

    ~GroundwaterFlowProcess()
    {
    }

private:
    ProcessVariable const* _hydraulic_head = nullptr;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
