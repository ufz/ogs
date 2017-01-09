/**
 * @file CsvInterface.h
 * @author Karsten Rink
 * @date 2015-03-25
 * @brief Definition of the CsvInterface class.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef CSVINTERFACE_H_
#define CSVINTERFACE_H_

#include <logog/include/logog.hpp>
#include <boost/any.hpp>

#include <array>
#include <fstream>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <typeinfo>
#include <vector>


#include "BaseLib/StringTools.h"
#include "BaseLib/IO/Writer.h"

namespace GeoLib {
    class Point;
}

namespace FileIO {

/**
 * Interface for reading CSV file formats.
 */
class CsvInterface  : public BaseLib::IO::Writer
{

public:
    /// Contructor (only needed for writing files)
    CsvInterface();

    /// Returns the number of vectors currently staged for writing.
    std::size_t getNArrays() const { return _vec_names.size(); }

    /// Adds an index vector of size s to the CSV file
    void addIndexVectorForWriting(std::size_t s);

    /// Stores if the CSV file to be written should include a header or not.
    void setCsvHeader(bool write_header) { _writeCsvHeader = write_header; }

    /// Adds a data vector to the CSV file. All data vectors have to have the same size.
    /// Vectors will be written in the same sequence they have been added to the interface.
    template<typename T>
    bool addVectorForWriting(std::string const& vec_name, std::vector<T> const& vec)
    {
        static_assert(std::is_same<T, std::string>::value
                || std::is_same<T, double>::value
                || std::is_same<T, int>::value,
                "CsvInterface can only write vectors of strings, doubles or ints.");

        if (!_data.empty())
        {
            std::size_t const vec_size (getVectorSize(0));
            if (vec_size != vec.size())
            {
                ERR ("Vector size does not match existing data (should be %d).", vec_size);
                return false;
            }
        }

        _vec_names.push_back(vec_name);
        _data.push_back(vec);
        return true;
    }

    /// Writes the CSV file.
    bool write();

    /**
     * Reads 3D points from a CSV file. It is assumed that the file has a header
     * specifying a name for each of the columns. The first three columns will be
     * interpreted as x-, y- and z-coordinate, respectively.
     * \param fname    Name of the file to be read
     * \param delim    Deliminator, default is ','
     * \param points   A vector containing the 3D points read from the file
     * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
     */
    static int readPoints(std::string const& fname, char delim,
                          std::vector<GeoLib::Point*> &points);

    /**
     * Reads 3D points from a CSV file. It is assumed that the file has a header
     * specifying a name for each of the columns. The columns specified in the
     * function call will be used for reading x-, y- and z-coordinates,
     * respectively If z_column_name is an empty string or not given at all, all
     * z-coordinates will be set to zero.
     * \param fname           Name of the file to be read
     * \param delim           Deliminator, default is ','
     * \param points          A vector containing the 3D points read from the file
     * \param x_column_name   Name of the column to be interpreted as x-coordinate
     * \param y_column_name   Name of the column to be interpreted as y-coordinate
     * \param z_column_name   Name of the column to be interpreted as z-coordinate
     * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
     */
    static int readPoints(std::string const& fname, char delim,
                          std::vector<GeoLib::Point*> &points,
                          std::string const& x_column_name,
                          std::string const& y_column_name,
                          std::string const& z_column_name = "");

    /**
     * Reads 3D points from a headerless CSV file, so columns for x-, y- and
     * z-coordinates have to be specified using indices (starting with 0).
     * If z_column_idx is not given (or set to numeric_limits::max()), all
     * z-coordinates will be set to zero.
     * \param fname          Name of the file to be read
     * \param delim          Deliminator, default is ','
     * \param points         A vector containing the 3D points read from the file
     * \param x_column_idx   Index of the column to be interpreted as x-coordinate
     * \param y_column_idx   Index of the column to be interpreted as y-coordinate
     * \param z_column_idx   Index of the column to be interpreted as z-coordinate
     * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
     */
    static int readPoints(std::string const& fname, char delim,
                          std::vector<GeoLib::Point*> &points,
                          std::size_t x_column_idx,
                          std::size_t y_column_idx,
                          std::size_t z_column_idx = std::numeric_limits<std::size_t>::max());

    /**
     * Reads a column of the given name from a CSV file.
     * \param fname        Name of the file to be read
     * \param delim        Deliminator, default is ','
     * \param data_array   A vector containing the data read from the file
     * \param column_name  The column's name to read
     * \return An error code (0 = ok, 0<i<max = number of skipped lines, -1 error reading file)
     */
    template <typename T>
    static int readColumn(std::string const& fname, char delim,
                          std::vector<T> &data_array,
                          std::string const& column_name)
    {
        std::ifstream in(fname.c_str());
        if (!in.is_open()) {
            ERR ("CsvInterface::readColumn(): Could not open file %s.", fname.c_str());
            return -1;
        }

        std::string line;
        getline(in, line);
        std::size_t const column_idx = CsvInterface::findColumn(line, delim, column_name);
        if (column_idx == std::numeric_limits<std::size_t>::max())
        {
            ERR ("Column \"%s\" not found in file header.", column_name.c_str());
            return -1;
        }
        return readColumn<T>(in, delim, data_array, column_idx);
    }

    template <typename T>
    static int readColumn(std::string const& fname, char delim,
                              std::vector<T> &data_array,
                              std::size_t column_idx)
    {
        std::ifstream in(fname.c_str());
        if (!in.is_open()) {
            ERR ("CsvInterface::readColumn(): Could not open file %s.", fname.c_str());
            return -1;
        }
        return readColumn<T>(in, delim, data_array, column_idx);
    }

private:
    /// Actual point reader for public methods
    static int readPoints(std::ifstream &in, char delim,
                          std::vector<GeoLib::Point*> &points,
                          std::array<std::size_t, 3> const& column_idx);

    /// Actual column reader for public methods
    template <typename T>
    static int readColumn(std::ifstream &in, char delim,
                       std::vector<T> &data_array,
                       std::size_t column_idx)
    {
        std::string line;
        std::size_t line_count(0);
        std::size_t error_count(0);
        std::list<std::string>::const_iterator it;
        while ( getline(in, line) )
        {
            line_count++;
            std::list<std::string> const fields = BaseLib::splitString(line, delim);

            if (fields.size() < column_idx+1)
            {
                ERR ("Line %d contains not enough columns of data. Skipping line...", line_count);
                error_count++;
                continue;
            }
            it = fields.begin();
            std::advance(it, column_idx);

            std::istringstream stream(*it);
            T value;
            if (!(stream >> value))
            {
                ERR("Error reading value in line %d.", line_count);
                error_count++;
                continue;
            }

            data_array.push_back(value);
        }
        return error_count;
    }

    /// Returns the number of the column with column_name (or std::numeric_limits::max() if no such column has been found).
    static std::size_t findColumn(std::string const& line, char delim, std::string const& column_name);

    /// Returns the size of the vector with the given index
    std::size_t getVectorSize(std::size_t idx) const;

    /**
     * Writes a value from a vector to the file.
     * \param vec_idx     Index of the vector
     * \param in_vec_idx  Entry in the selected vector
     */
    void writeValue(std::size_t vec_idx, std::size_t in_vec_idx);

    bool _writeCsvHeader;
    std::vector<std::string> _vec_names;
    std::vector< boost::any > _data;
};

} // FileIO

#endif /* CSVINTERFACE_H_ */
