/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-18
 * \brief  Definition of the XMLInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
class XMLInterface : public BaseLib::IO::Writer
{
public:
    XMLInterface();
    ~XMLInterface() override = default;

    void setNameForExport(std::string const& name) { _exportName = name; }
    virtual bool readFile(std::string const& fname) = 0;

protected:
    std::string _exportName;
};

}
}
