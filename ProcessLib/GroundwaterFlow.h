/**
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOW_H_
#define PROCESS_LIB_GROUNDWATERFLOW_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"
#include "OGS/ProcessVariable.h"

namespace ProcessLib
{

class GroundwaterFlowProcess : public Process
{
    using ConfigTree = boost::property_tree::ptree;

public:
    GroundwaterFlowProcess(MeshLib::Mesh const& mesh,
            std::vector<OGS::ProcessVariable> const& variables,
            ConfigTree const& config)
        : Process(mesh)
    {
        DBUG("Create GroundwaterFlowProcess.");

        // Find the corresponding process variable.
        std::string const name = config.get<std::string>("process_variable");

        auto const& variable = std::find_if(variables.cbegin(), variables.cend(),
                [&name](OGS::ProcessVariable const& v) {
                    return v.getName() == name;
                });

        if (variable == variables.end())
            ERR("Expected field \'%s\' not found in provided fields.",
                name.c_str());

        DBUG("Associate hydraulic_head with field named \'%s\'.",
            name.c_str());
        _hydraulic_head = &*variable;

    }

    void initialize()
    {
    }

private:
    OGS::ProcessVariable const* _hydraulic_head = nullptr;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_H_
