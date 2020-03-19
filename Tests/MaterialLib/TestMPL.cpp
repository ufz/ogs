/**
 * \file
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TestMPL.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/CreateMedium.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/Parameter.h"

std::unique_ptr<MPL::Medium> createTestMaterial(std::string const& xml)
{
    auto const ptree = readXml(xml.c_str());
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& config = conf.getConfigSubtree("medium");
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;

    return MPL::createMedium(config, parameters, nullptr, curves);
}
