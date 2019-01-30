/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FunctionParameter.h"

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "Function");

    std::vector<std::string> vec_expressions;

    //! \ogs_file_param{prj__parameters__parameter__Function__expression}
    for (auto p : config.getConfigSubtreeList("expression"))
    {
        std::string const expression_str = p.getValue<std::string>();
        vec_expressions.emplace_back(expression_str);
    }

    return std::make_unique<FunctionParameter<double>>(
        name, mesh, vec_expressions);
}

}  // ProcessLib
