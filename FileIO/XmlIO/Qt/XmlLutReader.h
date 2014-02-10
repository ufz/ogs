/**
 * \file
 * \author Karsten Rink
 * \date   2011-01-30
 * \brief  Definition of the XmlLutReader class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLLUTREADER_H
#define XMLLUTREADER_H

#include "VtkColorLookupTable.h"

#include <QFile>
#include <QtXml/QDomDocument>

// ThirdParty/logog
#include "logog/include/logog.hpp"


namespace FileIO
{

/**
 * \brief Reader for vtk-Lookup-Tables (in XML / ParaView Format)
 */
class XmlLutReader
{
public:
	static VtkColorLookupTable* readFromFile(const QString &fileName)
	{
		VtkColorLookupTable* lut = VtkColorLookupTable::New();

		QFile* file = new QFile(fileName);
		if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
		{
			ERR("XmlLutReader::readFromFile(): Can't open xml-file %s.", fileName.data());
			delete file;
			return NULL;
		}

		QDomDocument doc("ColorMap");
		doc.setContent(file);
		QDomElement docElement = doc.documentElement();
		if (docElement.nodeName().compare("ColorMap"))
		{
			ERR("XmlLutReader::readFromFile(): Unexpected XML root.");
			delete file;
			return NULL;
		}

		if (docElement.hasAttribute("interpolation"))
		{
			if (docElement.attribute("interpolation").compare("Linear") == 0)
				lut->setInterpolationType(VtkColorLookupTable::LUTType::LINEAR);
			else if (docElement.attribute("interpolation").compare("Exponential") == 0)
				lut->setInterpolationType(VtkColorLookupTable::LUTType::EXPONENTIAL);
			else
				lut->setInterpolationType(VtkColorLookupTable::LUTType::NONE);
		}
		else // default
			lut->setInterpolationType(VtkColorLookupTable::LUTType::NONE);

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

				unsigned char a[4] = { r, g, b, o };
				lut->setColor(value, a);
			}
			point = point.nextSiblingElement();
		}

		lut->SetTableRange(range[0], range[1]);

		delete file;

		return lut;
	};


};

 }

#endif // XMLLUTREADER_H
