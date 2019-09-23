/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateJacobianAssembler.h"
#include "BaseLib/Error.h"

#include "AnalyticalJacobianAssembler.h"
#include "CentralDifferencesJacobianAssembler.h"
#include "CompareJacobiansJacobianAssembler.h"

namespace ProcessLib
{
std::unique_ptr<AbstractJacobianAssembler> createJacobianAssembler(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return std::make_unique<AnalyticalJacobianAssembler>();
    }

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__type}
    auto const type = config->peekConfigParameter<std::string>("type");

    if (type == "Analytical") {
        config->ignoreConfigParameter("type");
        return std::make_unique<AnalyticalJacobianAssembler>();
    }
    if (type == "CentralDifferences")
    {
        return createCentralDifferencesJacobianAssembler(*config);
    }
    if (type == "CompareJacobians")
    {
        return createCompareJacobiansJacobianAssembler(*config);
    }

    OGS_FATAL("Unknown Jacobian assembler type: `%s'.", type.c_str());
}
}  // namespace ProcessLib
