/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Implementation of the XmlStnInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "XmlStnInterface.h"

#include <QFile>
#include <QtXml/QDomDocument>
#include <limits>
#include <rapidxml.hpp>

#include "BaseLib/DateTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/StationBorehole.h"

namespace GeoLib
{
namespace IO
{
XmlStnInterface::XmlStnInterface(GeoLib::GEOObjects& geo_objs)
    : XMLQtInterface("OpenGeoSysSTN.xsd"), _geo_objs(geo_objs)
{
}

int XmlStnInterface::readFile(const QString& fileName)
{
    if (XMLQtInterface::readFile(fileName) == 0)
    {
        return 0;
    }

    QDomDocument doc("OGS-STN-DOM");
    doc.setContent(getContent());
    QDomElement docElement =
        doc.documentElement();  // root element, used for identifying file-type
    if (docElement.nodeName().compare("OpenGeoSysSTN"))
    {
        ERR("XmlStnInterface::readFile(): Unexpected XML root.");
        return 0;
    }

    QDomNodeList lists = docElement.childNodes();
    for (int i = 0; i < lists.count(); i++)
    {
        // read all the station lists
        QDomNodeList stationList = lists.at(i).childNodes();
        std::vector<GeoLib::Point*> stations;
        std::string stnName("[NN]");

        for (int j = 0; j < stationList.count(); j++)
        {
            const QDomNode station_node(stationList.at(j));
            const QString station_type(station_node.nodeName());
            if (station_type.compare("name") == 0)
            {
                stnName = station_node.toElement().text().toStdString();
            }
            else if (station_type.compare("stations") == 0)
            {
                readStations(station_node, stations, fileName.toStdString());
            }
            else if (station_type.compare("boreholes") == 0)
            {
                readStations(station_node, stations, fileName.toStdString());
            }
        }

        if (!stations.empty())
        {
            _geo_objs.addStationVec(std::move(stations), stnName);
        }
    }

    return 1;
}

void XmlStnInterface::readStations(const QDomNode& stationsRoot,
                                   std::vector<GeoLib::Point*>& stations,
                                   const std::string& station_file_name)
{
    QDomElement station = stationsRoot.firstChildElement();
    while (!station.isNull())
    {
        if (station.hasAttribute("id") && station.hasAttribute("x") &&
            station.hasAttribute("y"))
        {
            std::string stationName("[NN]");
            std::string sensor_data_file_name;
            std::string boreholeDate("0000-00-00");
            double boreholeDepth(0.0);
            double stationValue(0.0);

            QDomNodeList stationFeatures = station.childNodes();
            for (int i = 0; i < stationFeatures.count(); i++)
            {
                // check for general station features
                const QDomNode feature_node(stationFeatures.at(i));
                const QString feature_name(feature_node.nodeName());
                const QString element_text(feature_node.toElement().text());
                if (feature_name.compare("name") == 0)
                {
                    stationName = element_text.toStdString();
                }
                if (feature_name.compare("sensordata") == 0)
                {
                    sensor_data_file_name = element_text.toStdString();
                    /* add other station features here */

                    // check for general borehole features
                }
                else if (feature_name.compare("value") == 0)
                {
                    stationValue = element_text.toDouble();
                }
                else if (feature_name.compare("bdepth") == 0)
                {
                    boreholeDepth = element_text.toDouble();
                }
                else if (feature_name.compare("bdate") == 0)
                {
                    boreholeDate = element_text.toStdString();
                }
                /* add other borehole features here */
            }

            double zVal = (station.hasAttribute("z"))
                              ? station.attribute("z").toDouble()
                              : 0.0;

            if (station.nodeName().compare("station") == 0)
            {
                GeoLib::Station* s =
                    new GeoLib::Station(station.attribute("x").toDouble(),
                                        station.attribute("y").toDouble(),
                                        zVal,
                                        stationName);
                s->setStationValue(stationValue);
                if (!sensor_data_file_name.empty())
                {
                    s->addSensorDataFromCSV(BaseLib::joinPaths(
                        station_file_name, sensor_data_file_name));
                }
                stations.push_back(s);
            }
            else if (station.nodeName().compare("borehole") == 0)
            {
                GeoLib::StationBorehole* s =
                    GeoLib::StationBorehole::createStation(
                        stationName,
                        station.attribute("x").toDouble(),
                        station.attribute("y").toDouble(),
                        zVal,
                        boreholeDepth,
                        boreholeDate);
                s->setStationValue(stationValue);
                /* add stratigraphy to the borehole */
                for (int j = 0; j < stationFeatures.count(); j++)
                {
                    if (stationFeatures.at(j).nodeName().compare("strat") == 0)
                    {
                        this->readStratigraphy(stationFeatures.at(j), s);
                    }
                }

                stations.push_back(s);
            }
        }
        else
        {
            WARN(
                "XmlStnInterface::readStations(): Attribute missing in "
                "<station> tag.");
        }
        station = station.nextSiblingElement();
    }
}

void XmlStnInterface::readStratigraphy(const QDomNode& stratRoot,
                                       GeoLib::StationBorehole* borehole)
{
    // borehole->addSoilLayer((*borehole)[0], (*borehole)[1], (*borehole)[2],
    // "");
    double depth_check((*borehole)[2]);
    QDomElement horizon = stratRoot.firstChildElement();
    while (!horizon.isNull())
    {
        if (horizon.hasAttribute("id") && horizon.hasAttribute("x") &&
            horizon.hasAttribute("y") && horizon.hasAttribute("z"))
        {
            std::string horizonName("[NN]");

            QDomNodeList horizonFeatures = horizon.childNodes();
            for (int i = 0; i < horizonFeatures.count(); i++)
            {
                if (horizonFeatures.at(i).nodeName().compare("name") == 0)
                {
                    horizonName =
                        horizonFeatures.at(i).toElement().text().toStdString();
                }
            }
            /* add other horizon features here */

            double depth(horizon.attribute("z").toDouble());
            if (std::abs(depth - depth_check) >
                std::numeric_limits<double>::
                    epsilon())  // skip soil-layer if its thickness is zero
            {
                borehole->addSoilLayer(horizon.attribute("x").toDouble(),
                                       horizon.attribute("y").toDouble(),
                                       depth,
                                       horizonName);
                depth_check = depth;
            }
            else
            {
                WARN(
                    "XmlStnInterface::readStratigraphy(): Skipped layer '{:s}' "
                    "in borehole '{:s}' because of thickness 0.0.",
                    horizonName, borehole->getName());
            }
        }
        else
        {
            WARN(
                "XmlStnInterface::readStratigraphy(): Attribute missing in "
                "<horizon> tag.");
        }
        horizon = horizon.nextSiblingElement();
    }
}

bool XmlStnInterface::write()
{
    if (export_name.empty())
    {
        ERR("XmlStnInterface::write(): No station list specified.");
        return false;
    }

    out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";  // xml
                                                                 // definition
    out << "<?xml-stylesheet type=\"text/xsl\" "
           "href=\"OpenGeoSysSTN.xsl\"?>\n\n";  // stylefile definition

    QDomDocument doc("OGS-STN-DOM");
    QDomElement root = doc.createElement("OpenGeoSysSTN");
    root.setAttribute("xmlns:ogs", "http://www.opengeosys.org");
    root.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");

    const std::vector<GeoLib::Point*>* stations(
        _geo_objs.getStationVec(export_name));
    bool const is_borehole =
        dynamic_cast<GeoLib::StationBorehole*>((*stations)[0]);

    doc.appendChild(root);
    QDomElement stationListTag = doc.createElement("stationlist");
    root.appendChild(stationListTag);

    QDomElement listNameTag = doc.createElement("name");
    stationListTag.appendChild(listNameTag);
    QDomText stationListNameText =
        doc.createTextNode(QString::fromStdString(export_name));
    listNameTag.appendChild(stationListNameText);
    QString listType = is_borehole ? "boreholes" : "stations";
    QDomElement stationsTag = doc.createElement(listType);
    stationListTag.appendChild(stationsTag);

    bool useStationValue(false);
    double sValue =
        static_cast<GeoLib::Station*>((*stations)[0])->getStationValue();
    std::size_t nStations(stations->size());
    for (std::size_t i = 1; i < nStations; i++)
    {
        if ((static_cast<GeoLib::Station*>((*stations)[i])->getStationValue() -
             sValue) < std::numeric_limits<double>::epsilon())
        {
            useStationValue = true;
            break;
        }
    }

    for (std::size_t i = 0; i < nStations; i++)
    {
        QString stationType = is_borehole ? "borehole" : "station";
        QDomElement stationTag = doc.createElement(stationType);
        stationTag.setAttribute("id", QString::number(i));
        stationTag.setAttribute(
            "x", QString::number((*(*stations)[i])[0], 'f',
                                 std::numeric_limits<double>::max_digits10));
        stationTag.setAttribute(
            "y", QString::number((*(*stations)[i])[1], 'f',
                                 std::numeric_limits<double>::max_digits10));
        stationTag.setAttribute(
            "z", QString::number((*(*stations)[i])[2], 'f',
                                 std::numeric_limits<double>::max_digits10));
        stationsTag.appendChild(stationTag);

        QDomElement stationNameTag = doc.createElement("name");
        stationTag.appendChild(stationNameTag);
        QDomText stationNameText = doc.createTextNode(QString::fromStdString(
            static_cast<GeoLib::Station*>((*stations)[i])->getName()));
        stationNameTag.appendChild(stationNameText);

        if (useStationValue)
        {
            QDomElement stationValueTag = doc.createElement("value");
            stationTag.appendChild(stationValueTag);
            QDomText stationValueText = doc.createTextNode(
                QString::number(static_cast<GeoLib::Station*>((*stations)[i])
                                    ->getStationValue()));
            stationValueTag.appendChild(stationValueText);
        }

        if (is_borehole)
        {
            writeBoreholeData(
                doc, stationTag,
                static_cast<GeoLib::StationBorehole*>((*stations)[i]));
        }
    }

    std::string xml = doc.toString().toStdString();
    out << xml;
    return true;
}

void XmlStnInterface::writeBoreholeData(QDomDocument& doc,
                                        QDomElement& boreholeTag,
                                        GeoLib::StationBorehole* borehole) const
{
    QDomElement stationDepthTag = doc.createElement("bdepth");
    boreholeTag.appendChild(stationDepthTag);
    QDomText stationDepthText =
        doc.createTextNode(QString::number(borehole->getDepth(), 'f'));
    stationDepthTag.appendChild(stationDepthText);
    if (std::abs(borehole->getDate()) > 0)
    {
        QDomElement stationDateTag = doc.createElement("bdate");
        boreholeTag.appendChild(stationDateTag);
        QDomText stationDateText = doc.createTextNode(
            QString::fromStdString(BaseLib::date2string(borehole->getDate())));
        stationDateTag.appendChild(stationDateText);
    }

    std::vector<GeoLib::Point*> profile = borehole->getProfile();
    std::vector<std::string> soilNames = borehole->getSoilNames();
    std::size_t nHorizons(profile.size());

    if (nHorizons > 1)
    {
        QDomElement stratTag = doc.createElement("strat");
        boreholeTag.appendChild(stratTag);

        for (std::size_t j = 1; j < nHorizons;
             j++)  /// the first entry in the profile vector is just the
                   /// position of the borehole
        {
            QDomElement horizonTag = doc.createElement("horizon");
            horizonTag.setAttribute("id", QString::number(j));
            horizonTag.setAttribute("x",
                                    QString::number((*profile[j])[0], 'f'));
            horizonTag.setAttribute("y",
                                    QString::number((*profile[j])[1], 'f'));
            horizonTag.setAttribute("z",
                                    QString::number((*profile[j])[2], 'f'));
            stratTag.appendChild(horizonTag);
            QDomElement horizonNameTag = doc.createElement("name");
            horizonTag.appendChild(horizonNameTag);
            QDomText horizonNameText =
                doc.createTextNode(QString::fromStdString(soilNames[j]));
            horizonNameTag.appendChild(horizonNameText);
        }
    }
}
}  // namespace IO
}  // namespace GeoLib
