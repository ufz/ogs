/**
 * \file
 * \author Lars Bilke
 * \date   2012-02-13
 * \brief  Definition of the Writer class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <filesystem.h>
#include <string>
#include <sstream>

namespace BaseLib
{
namespace IO
{

/// @brief Base class which enables writing an object to string, stringstream
/// or file. Also formatting (precision, scientific notation of decimal values)
/// can be set.
///
/// When subclassing you only need to implement void write() in which you have
/// to write to _out.
class Writer
{
public:
    Writer();
    virtual ~Writer() = default;

    /// @brief Writes the object to a string.
    std::string writeToString();

    /// @brief Writes the object to the given file.
    int writeToFile(std::filesystem::path const& file_path);

    /// @brief Sets the decimal precision.
    void setPrecision(unsigned int precision);

protected:
    /// @brief Writes the object to the internal stream.
    /// This method must be implemented by a subclass.
    /// The implementation should return true on success, else false
    virtual bool write() = 0;

    /// @brief The stream to write to.
    std::stringstream _out;

private:

};

} // namespace IO
} // namespace BaseLib
