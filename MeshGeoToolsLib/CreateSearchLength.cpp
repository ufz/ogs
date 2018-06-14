/**
 * @file
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "CreateSearchLength.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MeshGeoToolsLib/HeuristicSearchLength.h"
#include "MeshGeoToolsLib/SearchLength.h"

namespace MeshGeoToolsLib
{
std::unique_ptr<MeshGeoToolsLib::SearchLength> createSearchLengthAlgorithm(
    BaseLib::ConfigTree const& external_config, MeshLib::Mesh const& mesh)
{
    boost::optional<BaseLib::ConfigTree> config =
        //! \ogs_file_param{prj__search_length_algorithm}
        external_config.getConfigSubtreeOptional("search_length_algorithm");

    if (!config)
    {
        return std::make_unique<MeshGeoToolsLib::SearchLength>();
    }

    //! \ogs_file_param{prj__search_length_algorithm__type}
    std::string const type = config->getConfigParameter<std::string>("type");

    //! \ogs_file_param_special{prj__search_length_algorithm__fixed}
    if (type == "fixed")
    {
        //! \ogs_file_param{prj__search_length_algorithm__fixed__value}
        double const length = config->getConfigParameter<double>("value");
        return std::make_unique<MeshGeoToolsLib::SearchLength>(length);
    }
    if (type == "heuristic")
    {
        //! \ogs_file_param_special{prj__search_length_algorithm__heuristic}
        return std::make_unique<HeuristicSearchLength>(mesh);
    }
    OGS_FATAL("Unknown search length algorithm type '%s'.", type.c_str());
}

}  // end namespace MeshGeoToolsLib
