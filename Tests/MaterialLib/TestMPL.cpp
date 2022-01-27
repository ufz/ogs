/**
 * \file
 * \author Norbert Grunwald
 * \date   Oct 22, 2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TestMPL.h"

#include <sstream>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/CreateMedium.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "Tests/TestTools.h"

namespace Tests
{
std::unique_ptr<MPL::Medium> createTestMaterial(
    std::string const& xml, int const geometry_dimension,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    auto ptree = Tests::readXml(xml.c_str());
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& config = conf.getConfigSubtree("medium");
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>>
        curves;

    return MPL::createMedium(0 /*material id*/, geometry_dimension, config,
                             parameters, local_coordinate_system, curves);
}

std::unique_ptr<MaterialPropertyLib::Property> createTestProperty(
    const char xml[],
    std::function<std::unique_ptr<MaterialPropertyLib::Property>(
        BaseLib::ConfigTree const& config)>
        createProperty)
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("property");
    // Parsing the property name:
    auto const property_name =
        sub_config.getConfigParameter<std::string>("name");

    return createProperty(sub_config);
}

std::string makeConstantPropertyElement(std::string const name,
                                        double const value)
{
    std::stringstream ss;

    ss << "<property>\n";
    ss << "    <name>" << name << "</name>\n";
    ss << "    <type>Constant</type>\n";
    ss << "    <value> " << value << " </value>\n";
    ss << "</property>\n";

    return ss.str();
}
}  // namespace Tests
