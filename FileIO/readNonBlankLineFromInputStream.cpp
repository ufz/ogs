/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file readNonBlankLineFromInputStream.cpp
 *
 *  Created on 2011-04-19 by Thomas Fischer
 */

#include "readNonBlankLineFromInputStream.h"

std::string readNonBlankLineFromInputStream(std::istream & in)
{
	std::string line;

	bool not_finished (true);
	while (not_finished)
	{
		// read line
		getline(in, line);
		if (!in.fail())
		{
			// skip initial space characters
			std::string::size_type i (line.find_first_not_of(" ", 0));
			// search comment symbol ;
			std::string::size_type j (line.find(";", i));
			if (j == i) // first non space character is equal to the comment symbol
				not_finished = true;
			else
			{
				if ((i != std::string::npos))
					// cut string from first non blank space character until the first comment character
					line = line.substr(i, j - i);

				// remove last blank spaces
				i = line.find_last_not_of(" ");
				if (i != std::string::npos)
					line = line.substr(0, i + 1);

				not_finished = false;
			}
		}
		else
		{
			line = "";
			not_finished = false;
		}
	}
	return line;
}
