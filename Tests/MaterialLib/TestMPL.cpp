/**
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TestMPL.h"

MPL::Medium createTestMaterial(std::string const& xml)
{
    auto const ptree = readXml(xml.c_str());
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& config = conf.getConfigSubtree("medium");

    return MPL::Medium(config);
}
