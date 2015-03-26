/**
 * @file CsvInterface.h
 * @author Karsten Rink
 * @date 2015-03-25
 * @brief Definition of the CsvInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef CSVINTERFACE_H_
#define CSVINTERFACE_H_

#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <limits>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "StringTools.h"

namespace GeoLib {
	class Point;
}

namespace FileIO {

/**
 * Interface for reading CSV file formats.
 */
class CsvInterface {

public:
	/** 
	 * Reads 3D points from a CSV file. It is assumed that the file has a header
	 * specifying a name for each of the rows. The first three rows will be 
	 * interpreted as x-, y- and z-coordinate, respectively.
	 * \param fname    Name of the file to be read
	 * \param delim    Deliminator, default is ','
	 * \param points   A vector containing the 3D points read from the file
	 * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
	 */
	static int readPointsFromCSV(std::string const& fname, char delim,
	                             std::vector<GeoLib::Point*> &points);

	/** 
	 * Reads 3D points from a CSV file. It is assumed that the file has a header
	 * specifying a name for each of the rows. The rows specified in the function
	 * call will be used for reading x-, y- and z-coordinates, respectively
	 * If z_row_name is an empty string or not given at all, all z-coordinates
	 * will be set to zero.
	 * \param fname        Name of the file to be read
	 * \param delim        Deliminator, default is ','
	 * \param points       A vector containing the 3D points read from the file
	 * \param x_row_name   Name of the row to be interpreted as x-coordinate
	 * \param y_row_name   Name of the row to be interpreted as y-coordinate
	 * \param z_row_name   Name of the row to be interpreted as z-coordinate
	 * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
	 */
	static int readPointsFromCSV(std::string const& fname, char delim,
	                             std::vector<GeoLib::Point*> &points,
	                             std::string const& x_row_name,
	                             std::string const& y_row_name,
	                             std::string const& z_row_name);

	/**
	 * Reads a row of the given name from a CSV file.
	 * \param fname        Name of the file to be read
	 * \param delim        Deliminator, default is ','
	 * \param data_arary   A vector containing the data read from the file
	 * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
	 */
	template <typename T>
	static int readRowFromCSV(std::string const& fname, char delim,
	                          std::vector<T> &data_array,
	                          std::string const& row_name)
	{
		std::ifstream in(fname.c_str());
		if (!in.is_open()) {
			ERR ("CsvInterface::readPointsFromCSV(): Could not open file %s.", fname.c_str());
			return -1;
		}

		std::string line;
		getline(in, line);
		std::size_t const row_idx = CsvInterface::findRow(line, delim, row_name);
		if (row_idx == std::numeric_limits<std::size_t>::max())
		{
			ERR ("Row \"%s\" not found in file header.", row_name.c_str());
			return -1;
		}

		std::size_t line_count(0);
		std::size_t error_count(0);
		std::list<std::string>::const_iterator it;
		while ( getline(in, line) )
		{
			line_count++;
			std::list<std::string> const fields = BaseLib::splitString(line, delim);

			if (fields.size() < row_idx+1)
			{
				ERR ("Line %d contains not enough rows of data. Skipping line...", line_count);
				error_count++;
				continue;
			}
			it = fields.begin();
			std::advance(it, row_idx);

			std::istringstream stream(*it);
			T value;
			stream >> value;

			data_array.push_back(value);
		}
		return error_count;
	}

private:
	/// Returns the number of the row with row_name (or std::numeric_limits::max() if no such row has been found).
	static std::size_t findRow(std::string const& line, char delim, std::string const& row_name);
};

}

#endif /* CSVINTERFACE_H_ */
