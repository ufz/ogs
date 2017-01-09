/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlGspInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLGSPINTERFACE_H
#define XMLGSPINTERFACE_H

#include <vector>
#include <string>

#include <QString>

#include "BaseLib/IO/XmlIO/XMLInterface.h"
#include "BaseLib/IO/XmlIO/Qt/XMLQtInterface.h"

#include "Applications/DataHolderLib/Project.h"

namespace FileIO
{

/**
 * \brief Reads and writes project information to and from XML files.
 */
class XmlGspInterface : public BaseLib::IO::XMLInterface,
                        public BaseLib::IO::XMLQtInterface
{
public:
    XmlGspInterface(DataHolderLib::Project &project);

    virtual ~XmlGspInterface() {}

    /// Reads an xml-file containing a GeoSys project.
    /// Project files currently cover only geo-, msh- and station-data. This will be expanded in the future.
    int readFile(const QString &fileName);

    bool readFile(std::string const& fname) { return readFile(QString(fname.c_str())) != 0; }

    int writeToFile(const std::string& filename);

protected:
    bool write();

private:
    std::string _filename;
    DataHolderLib::Project& _project;
};

} // end namespace FileIO

#endif // XMLGSPINTERFACE_H
