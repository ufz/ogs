/**
 * @brief Functionality to build different search length algorithm objects from
 * given config.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include "BaseLib/ConfigTree.h"

#include "MeshGeoToolsLib/HeuristicSearchLength.h"
#include "MeshGeoToolsLib/SearchLength.h"

namespace MeshGeoToolsLib
{
MeshGeoToolsLib::SearchLength createSearchLengthAlgorithm(
    BaseLib::ConfigTree const& external_config, MeshLib::Mesh const& mesh)
{
    boost::optional<BaseLib::ConfigTree> config =
        //! \ogs_file_param{prj__search_length_algorithm}
        external_config.getConfigSubtreeOptional("search_length_algorithm");

    if (!config)
        return MeshGeoToolsLib::SearchLength();

    //! \ogs_file_param{prj__search_length_algorithm__type}
    std::string const type = config->getConfigParameter<std::string>("type");

    if (type == "fixed")
    {
        //! \ogs_file_param{prj__search_length_algorithm__value}
        double const length = config->getConfigParameter<double>("value");
        return MeshGeoToolsLib::SearchLength(length);
    }
    if (type == "heuristic")
    {
        return MeshGeoToolsLib::HeuristicSearchLength(mesh);
    }
    OGS_FATAL("Unknown search length algorithm type '%s'.", type.c_str());
}

}  // end namespace MeshGeoToolsLib
