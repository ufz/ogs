/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Color.cpp
 *
 * Created on 2010-06-16 by Karsten Rink
 */

#include <iostream>
#include <sstream>

#include "Color.h"
#include "StringTools.h"

namespace GeoLib {

Color* getRandomColor()
{
	return new Color((rand()%5)*50, (rand()%5)*50, (rand()%5)*50);
}

int readColorLookupTable(std::map<std::string, Color*> &colors, const std::string &filename)
{
	std::string id = "", line = "";

	std::ifstream in( filename.c_str() );

	if (!in.is_open())
	{
		std::cout << "Color::readLookupTable() - Could not open file..."  << std::endl;
		return 0;
	}

	while ( getline(in, line) )
	{
		std::list<std::string> fields = splitString(line, '\t');
		Color *c = new Color();

		if (fields.size()>=4)
		{
			id = fields.front();
			fields.pop_front();
			(*c)[0] = atoi(fields.front().c_str());
			fields.pop_front();
			(*c)[1] = atoi(fields.front().c_str());
			fields.pop_front();
			(*c)[2] = atoi(fields.front().c_str());
			colors.insert(std::pair<std::string, Color*>(id, c));
		}
		else delete c;
	}

	return 1;
}

const Color* getColor(const std::string &id, std::map<std::string, Color*> &colors)
{
	for (std::map<std::string, Color*>::const_iterator it=colors.begin(); it !=colors.end(); ++it)
	{
		if (id.compare(it->first) == 0)
			return it->second;
	}
	std::cout << "Key \"" << id << "\" not found in color lookup table..." << std::endl;
	Color* c = getRandomColor();
	colors.insert(std::pair<std::string, Color*>(id, c));
	return c;
}

const Color* getColor(double id, std::map<std::string, GeoLib::Color*> &colors)
{
	std::ostringstream stream;
	stream << id;
	return getColor(stream.str(), colors);
}

}
