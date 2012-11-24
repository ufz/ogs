/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file OptionsXMLReader.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "Options.h"

namespace BaseLib
{

/**
 * read a XML file and add contents to BaseLib::Options
 *
 * \param xml_file  a path to the XML file
 * \param propertoes BaseLib::Options object
 */
bool addXMLtoOptions(const std::string &xml_file, BaseLib::Options &properties);

}
