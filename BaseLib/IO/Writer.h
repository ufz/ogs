// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <sstream>
#include <string>
#include <string_view>

namespace BaseLib
{
namespace IO
{
/// @brief Base class which enables writing an object to string, stringstream
/// or file.
///
/// When subclassing you only need to implement void write() in which you have
/// to write to out.
class Writer
{
public:
    Writer();
    virtual ~Writer() = default;

    /// @brief Writes the object to a string.
    std::string writeToString();

protected:
    /// @brief Writes the object to the internal stream.
    /// This method must be implemented by a subclass.
    /// The implementation should return true on success, else false
    virtual bool write() = 0;

    /// @brief The stream to write to.
    std::ostringstream out;
};

/// \returns 0 if string is empty, or if there is an error, and 1 otherwise.
int writeStringToFile(std::string_view content,
                      std::filesystem::path const& file_path);
}  // namespace IO
}  // namespace BaseLib
