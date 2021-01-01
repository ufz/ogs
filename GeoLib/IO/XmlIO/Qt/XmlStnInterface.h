/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlStnInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/IO/XmlIO/XMLInterface.h"
#include "BaseLib/IO/XmlIO/Qt/XMLQtInterface.h"

namespace GeoLib {
class GEOObjects;
class Point;
class StationBorehole;

namespace IO
{

/**
 * \brief Reads and writes Observation Sites to and from XML files.
 */
class XmlStnInterface final : public BaseLib::IO::XMLInterface,
                              public BaseLib::IO::XMLQtInterface
{
public:
    explicit XmlStnInterface(GeoLib::GEOObjects& geo_objs);

    /// Reads an xml-file containing station object definitions into the GEOObjects used in the contructor (requires Qt)
    int readFile(const QString& fileName) override;

    bool readFile(std::string const& fname) override
    {
        return readFile(QString(fname.c_str())) != 0;
    }

protected:
    bool write() override;

private:
    /// Reads GeoLib::Station- or StationBorehole-objects from an xml-file
    void readStations  ( const QDomNode &stationsRoot, std::vector<GeoLib::Point*>* stations, const std::string &station_file_name);

    /// Writes borehole-specific data to a station-xml-file.
    void writeBoreholeData(QDomDocument &doc,
                           QDomElement &boreholeTag,
                           GeoLib::StationBorehole* borehole) const;

    /// Reads the stratigraphy of a borehole from an xml-file
    void readStratigraphy( const QDomNode &stratRoot, GeoLib::StationBorehole*  borehole );

    GeoLib::GEOObjects& _geo_objs;
};

} // end namespace IO
} // end namespace GeoLib
