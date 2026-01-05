// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

#include "BaseLib/IO/Writer.h"

namespace BaseLib
{
namespace IO
{
/**
 * \brief Base class for writing any information to and from XML files.
 */
struct XMLInterface : public BaseLib::IO::Writer
{
    virtual bool readFile(std::string const& fname) = 0;

    std::string export_name = {};
};

}  // namespace IO
}  // namespace BaseLib
