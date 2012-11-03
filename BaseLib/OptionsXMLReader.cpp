/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file OptionsXMLReader.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include "OptionsXMLReader.h"

#include <fstream>
#include <map>
#include <vector>

#include "logog.hpp"
#include "RapidXML/rapidxml.hpp"

#include "StringTools.h"

namespace BaseLib
{

void addXMLtoOptions(rapidxml::xml_node<>* e_root, BaseLib::OptionGroup &properties)
{
    BaseLib::OptionGroup* optRoot = properties.addSubGroup(e_root->name());

    // attributes
    for (const rapidxml::xml_attribute<>* att = e_root->first_attribute(); att != 0; att = att->next_attribute())
    {
        optRoot->addOption(att->name(), att->value());
    }

    // element
    for (rapidxml::xml_node<>* e=e_root->first_node(); e!=0; e=e->next_sibling())
    {
        if (e->first_node()!=0 || e->first_attribute()!=0) {
            addXMLtoOptions(e, *optRoot);
        } else if (e->value() != 0){
            optRoot->addOption(e->name(), e->value());
        }
    }
}

bool addXMLtoOptions(const std::string &file_name, BaseLib::Options &properties)
{
    // 
	std::ifstream in(file_name.c_str());
	if (in.fail())
	{
		std::cout << "\nBaseLib::addXMLtoOptions() - Can't open xml-file." << std::endl;
		return false;
	}

	in.seekg(0, std::ios::end);
    const std::size_t length = in.tellg();
	in.seekg(0, std::ios::beg);
	char* buffer = new char[length+1];
	in.read(buffer, length);
	buffer[in.gcount()] = '\0';
	in.close();

    rapidxml::xml_document<> doc;
    try {
        doc.parse<0>(buffer);
    } catch (rapidxml::parse_error &e) {
        LOGOG_CERR << "Error in reading a XML file " << file_name << " with error " << e.what() << std::endl;
        return false;
    }

    // read data
    addXMLtoOptions(doc.first_node(), properties);


    return true;
}

} //end
