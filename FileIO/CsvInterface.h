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

#include <vector>

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
	 * \param fname    Name of the file to be read.
	 * \param delim    Deliminator, default is ','
	 * \result A vector containing 3D points.
	 */
	static std::vector<GeoLib::Point*> readPointsFromCSV(std::string const& fname, char delim);
	
	/** 
	 * Reads 3D points from a CSV file. It is assumed that the file has a header
	 * specifying a name for each of the rows. The rows specified in the function 
	 * call will be used for reading x-, y- and z-coordinates, respectively
	 * \param fname      Name of the file to be read.
	 * \param delim      Deliminator, default is ','
	 * \param x_row_name Name of the row to be interpreted as x-coordinate
	 * \param y_row_name Name of the row to be interpreted as y-coordinate
	 * \param z_row_name Name of the row to be interpreted as z-coordinate
	 * \result A vector containing 3D points.
	 */
	static std::vector<GeoLib::Point*> readPointsFromCSV(std::string const& fname, char delim,
	                                                     std::string const& x_row_name,
	                                                     std::string const& y_row_name,
	                                                     std::string const& z_row_name);

private:
	/// Returns the number of the row with row_name (or std::numeric_limits::max() if no such row has been found).
	static std::size_t findRow(std::string const& line, char delim, std::string const& row_name);
};

}

#endif /* CSVINTERFACE_H_ */
