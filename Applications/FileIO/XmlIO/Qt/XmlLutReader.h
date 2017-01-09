/**
 * \file
 * \author Karsten Rink
 * \date   2011-01-30
 * \brief  Definition of the XmlLutReader class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLLUTREADER_H
#define XMLLUTREADER_H

#include <QFile>
#include <QtXml/QDomDocument>

#include <logog/include/logog.hpp>

#include "Applications/DataHolderLib/ColorLookupTable.h"


namespace FileIO
{

/**
 * \brief Reader for vtk-Lookup-Tables (in XML / ParaView Format)
 */
class XmlLutReader
{
public:
    static bool readFromFile(const QString &fileName, DataHolderLib::ColorLookupTable &lut)
    {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            ERR("XmlLutReader::readFromFile(): Can't open xml-file %s.", fileName.data());
            return false;
        }

        QDomDocument doc("ColorMap");
        doc.setContent(&file);
        QDomElement docElement = doc.documentElement();
        if (docElement.nodeName().compare("ColorMap"))
        {
            ERR("XmlLutReader::readFromFile(): Unexpected XML root.");
            file.close();
            return false;
        }

        if (docElement.hasAttribute("interpolation"))
        {
            if (docElement.attribute("interpolation").compare("Linear") == 0)
                lut.setInterpolationType(DataHolderLib::LUTType::LINEAR);
            else if (docElement.attribute("interpolation").compare("Exponential") == 0)
                lut.setInterpolationType(DataHolderLib::LUTType::EXPONENTIAL);
            else
                lut.setInterpolationType(DataHolderLib::LUTType::NONE);
        }
        else // default
            lut.setInterpolationType(DataHolderLib::LUTType::NONE);

        QDomElement point = docElement.firstChildElement();
        double range[2] = { point.attribute("x").toDouble(), point.attribute("x").toDouble() };

        while (!point.isNull())
        {
            if ((point.nodeName().compare("Point") == 0 )
                && point.hasAttribute("x")
                && point.hasAttribute("r")
                && point.hasAttribute("g")
                && point.hasAttribute("b"))
            {
                double value = point.attribute("x").toDouble();
                unsigned char r = static_cast<int>(255 * point.attribute("r").toDouble());
                unsigned char g = static_cast<int>(255 * point.attribute("g").toDouble());
                unsigned char b = static_cast<int>(255 * point.attribute("b").toDouble());
                unsigned char o = static_cast<int>(255 * (point.hasAttribute("o") ? point.attribute("o").toDouble() : 1));

                if (value < range[0])
                    range[0] = value;
                if (value > range[1])
                    range[1] = value;

                DataHolderLib::Color color = DataHolderLib::createColor(r, g, b, o);
                lut.setColor(value, color);
            }
            point = point.nextSiblingElement();
        }

        lut.setTableRange(range[0], range[1]);

        file.close();
        return true;
    };


};

 }

#endif // XMLLUTREADER_H
