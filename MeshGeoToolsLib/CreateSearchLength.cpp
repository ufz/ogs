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
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__search_length_algorithm}
        external_config.getConfigSubtreeOptional("search_length_algorithm");

    if (!config)
        return std::unique_ptr<MeshGeoToolsLib::SearchLength>{
            new MeshGeoToolsLib::SearchLength()};

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__search_length_algorithm__type}
    std::string const type = config->getConfigParameter<std::string>("type");

    if (type == "fixed")
    {
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__search_length_algorithm__fixed__value}
        double const length = config->getConfigParameter<double>("value");
        return std::unique_ptr<MeshGeoToolsLib::SearchLength>{
            new MeshGeoToolsLib::SearchLength(length)};
    }
    if (type == "heuristic")
    {
        return std::unique_ptr<MeshGeoToolsLib::HeuristicSearchLength>{
            new MeshGeoToolsLib::HeuristicSearchLength(mesh)};
    }
    OGS_FATAL("Unknown search length algorithm type '%s'.", type.c_str());
}

}  // end namespace MeshGeoToolsLib
