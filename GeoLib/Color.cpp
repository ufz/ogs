/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-16
 * \brief  Implementation of the Color class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <sstream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "Color.h"
#include "StringTools.h"

namespace GeoLib {

Color* getRandomColor()
{
	std::array<unsigned char, 3> col;
	col[0] = static_cast<unsigned char>((rand()%5)*50);
	col[1] = static_cast<unsigned char>((rand()%5)*50);
	col[2] = static_cast<unsigned char>((rand()%5)*50);
	return new Color(col);
}

int readColorLookupTable(std::map<std::string, Color*> &colors, const std::string &filename)
{
	std::string id = "", line = "";

	std::ifstream in( filename.c_str() );

	if (!in.is_open())
	{
		WARN("Color::readLookupTable() - Could not open file %s.", filename.c_str());
		return 0;
	}

	while ( getline(in, line) )
	{
		std::list<std::string> fields = BaseLib::splitString(line, '\t');

		if (fields.size()>=4)
		{
			Color *c = new Color();
			id = fields.front();
			fields.pop_front();
			(*c)[0] = atoi(fields.front().c_str());
			fields.pop_front();
			(*c)[1] = atoi(fields.front().c_str());
			fields.pop_front();
			(*c)[2] = atoi(fields.front().c_str());
			colors.insert(std::pair<std::string, Color*>(id, c));
		}
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
	WARN("Key \"%s\" not found in color lookup table.", id.c_str());
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
