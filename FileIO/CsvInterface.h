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

#include <array>
#include <fstream>
#include <limits>
#include <list>
#include <string>
#include <vector>

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
	                             std::string const& z_row_name = "");

	/** 
	 * Reads 3D points from a headerless CSV file, so rows for x-, y- and 
	 * z-coordinates have to be specified using indices (starting with 0).
	 * If z_row_idx is not given (or set to numeric_limits::max()), all 
	 * z-coordinates will be set to zero.
	 * \param fname       Name of the file to be read
	 * \param delim       Deliminator, default is ','
	 * \param points      A vector containing the 3D points read from the file
	 * \param x_row_idx   Index of the row to be interpreted as x-coordinate
	 * \param y_row_idx   Index of the row to be interpreted as y-coordinate
	 * \param z_row_idx   Index of the row to be interpreted as z-coordinate
	 * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
	 */
	static int readPointsFromCSV(std::string const& fname, char delim,
	                             std::vector<GeoLib::Point*> &points,
	                             std::size_t x_row_idx,
	                             std::size_t y_row_idx,
	                             std::size_t z_row_idx = std::numeric_limits<std::size_t>::max());

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
		return readRow<T>(in, delim, data_array, row_idx);
	}

	template <typename T>
	static int readRowFromCSV(std::string const& fname, char delim,
	                          std::vector<T> &data_array,
	                          std::size_t row_idx)
	{
		std::ifstream in(fname.c_str());
		if (!in.is_open()) {
			ERR ("CsvInterface::readPointsFromCSV(): Could not open file %s.", fname.c_str());
			return -1;
		}
		return readRow<T>(in, delim, data_array, row_idx);
	}

private:
	/// Actual point reader for public methods
	static int readPoints(std::ifstream &in, char delim,
	                      std::vector<GeoLib::Point*> &points,
	                      std::array<std::size_t, 3> const& row_idx);

	/// Actual row reader for public methods
	template <typename T>
	static int readRow(std::ifstream &in, char delim,
	                   std::vector<T> &data_array,
	                   std::size_t row_idx)
	{
		std::string line;
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

	/// Returns the number of the row with row_name (or std::numeric_limits::max() if no such row has been found).
	static std::size_t findRow(std::string const& line, char delim, std::string const& row_name);
};

}

#endif /* CSVINTERFACE_H_ */
