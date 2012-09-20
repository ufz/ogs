/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Writer.h
 *
 * Created on 2012-02-13 by Lars Bilke
 */

#ifndef WRITER_H
#define WRITER_H

#include <string>
#include <iostream>
#include <sstream>

namespace FileIO
{

/// @brief Base class which enables writing an object to string, stringstream
/// or file. Also formatting (precision, scientific notation of decimal values)
/// can be set.
///
/// When subclassing you only need to implement void write(std::ostream& stream).
class Writer
{
public:
	Writer();
	virtual ~Writer() {};

	/// @brief Writes the object to a string.
	std::string writeToString();

	/// @brief Writes the object to the given file.
	int writeToFile(std::string const& filename);

	/// @brief Sets the decimal precision.
	void setPrecision(unsigned int precision);

	/// @brief Sets the format (either ios::scientific or ios::fixed);
	void setFormat(std::ios_base::fmtflags flags);

protected:
	/// @brief Writes the object to the given stream.
	/// This method must be implemented by a subclass.
	virtual int write(std::ostream& stream) = 0;

	/// @brief The stream to write to.
	std::stringstream _out;

private:

};

} // namespace FileIO

#endif // WRITER_H
