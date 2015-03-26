/**
 * @file CsvInterface.cpp
 * @author Karsten Rink
 * @date 2015-03-25
 * @brief Implementation of the CsvInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CsvInterface.h"

#include <array>

#include "GeoLib/Point.h"

namespace FileIO {

int CsvInterface::readPointsFromCSV(std::string const& fname, char delim, 
                                    std::vector<GeoLib::Point*> &points)
{
	std::ifstream in(fname.c_str());

	if (!in.is_open()) {
		ERR ("CsvInterface::readPointsFromCSV(): Could not open file %s.", fname.c_str());
		return -1;
	}

	std::string line;
	getline(in, line);

	std::size_t line_count(0);
	std::size_t error_count(0);
	std::list<std::string>::const_iterator it;
	while ( getline(in, line) )
	{
		line_count++;
		std::list<std::string> const fields = BaseLib::splitString(line, delim);

		if (fields.size() < 3)
		{
			ERR ("Line %d contains not enough rows of data. Skipping line...", line_count);
			error_count++;
			continue;
		}
		it = fields.begin();
		std::array<double, 3> point;
		point[0] = strtod(it->c_str(), 0);
		point[1] = strtod((++it)->c_str(), 0);
		point[2] = strtod((++it)->c_str(), 0);
		points.push_back(new GeoLib::Point(point[0], point[1], point[2]));
	}
	return error_count;
}

int CsvInterface::readPointsFromCSV(std::string const& fname, char delim,
                                    std::vector<GeoLib::Point*> &points,
                                    std::string const& x_row_name,
                                    std::string const& y_row_name,
                                    std::string const& z_row_name = "")
{
	std::ifstream in(fname.c_str());
	std::array<std::string, 3> const row_names = { x_row_name, y_row_name, z_row_name };

	if (!in.is_open()) {
		ERR ("CsvInterface::readPointsFromCSV(): Could not open file %s.", fname.c_str());
		return -1;
	}

	std::string line;
	getline(in, line);
	std::array<std::size_t, 3> const row_idx = 
		{ CsvInterface::findRow(line, delim, x_row_name),
		  CsvInterface::findRow(line, delim, y_row_name),
		  (z_row_name.empty()) ? CsvInterface::findRow(line, delim, y_row_name) : CsvInterface::findRow(line, delim, z_row_name) };
	
	for (std::size_t i=0; i<3; ++i)
		if (row_idx[i] == std::numeric_limits<std::size_t>::max())
		{
			ERR ("Row \"%s\" not found in file header.", row_names[i].c_str());
			if (!(i == 2 || row_names[2].empty()))
				return -1;
		}

	std::array<std::size_t, 3> order = { 0, 1, 2 };
	std::sort(order.begin(), order.end(), 
		[&row_idx](std::size_t idx1, std::size_t idx2) {return row_idx[idx1] < row_idx[idx2];});
	std::array<std::size_t, 3> const row_advance = 
		{ row_idx[order[0]], 
		  row_idx[order[1]]-row_idx[order[0]], 
		  row_idx[order[2]]-row_idx[order[1]] };

	std::size_t line_count(0);
	std::size_t error_count(0);
	std::list<std::string>::const_iterator it;
	while ( getline(in, line) )
	{
		line_count++;
		std::list<std::string> const fields = BaseLib::splitString(line, delim);

		if (fields.size() < row_idx[order[2]]+1)
		{
			ERR ("Line %d contains not enough rows of data. Skipping line...", line_count);
			error_count++;
			continue;
		}

		std::array<double, 3> point;
		it = fields.begin();
		std::advance(it, row_advance[0]);
		point[0] = strtod(it->c_str(), 0);
		std::advance(it, row_advance[1]);
		point[1] = strtod(it->c_str(), 0);
		std::advance(it, row_advance[2]);
		point[2] = (z_row_name.empty()) ? 0 : strtod(it->c_str(), 0);
		points.push_back(new GeoLib::Point(point[0], point[1], point[2]));
	}
	return error_count;
}

std::size_t CsvInterface::findRow(std::string const& line, char delim, std::string const& row_name)
{
	std::list<std::string> const fields = BaseLib::splitString(line, delim);
	if (fields.size() < 3)
	{
		ERR (" The csv-file needs to contain at least three rows of data");
		return std::numeric_limits<std::size_t>::max();
	}

	std::size_t count(0);
	for (auto it = fields.cbegin(); it != fields.cend(); ++it)
	{
		if ((*it).compare(row_name) == 0)
			break;
		else
			count++;
	}

	if (count == fields.size())
		return std::numeric_limits<std::size_t>::max();

	return count;
}

} // end namespace FileIO
