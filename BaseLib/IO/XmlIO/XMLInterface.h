/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-18
 * \brief  Definition of the XMLInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
