/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file XmlLutReader.h
 *
 * Created on 2011-01-30 by Karsten Rink
 */

#ifndef XMLLUTREADER_H
#define XMLLUTREADER_H

#include "VtkColorLookupTable.h"

#include <QFile>
#include <QtXml/QDomDocument>

#include <iostream>

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
			std::cout << "XmlLutReader::readFromFile() - Can't open xml-file." << std::endl;
			delete file;
			return NULL;
		}

		QDomDocument doc("ColorMap");
		doc.setContent(file);
		QDomElement docElement = doc.documentElement();
		if (docElement.nodeName().compare("ColorMap"))
		{
			std::cout << "XmlLutReader::readFromFile() - Unexpected XML root." << std::endl;
			delete file;
			return NULL;
		}

		if (docElement.hasAttribute("interpolation"))
		{
			if (docElement.attribute("interpolation").compare("Linear") == 0)
				lut->setInterpolationType(VtkColorLookupTable::LINEAR);
			else if (docElement.attribute("interpolation").compare("Exponential") == 0)
				lut->setInterpolationType(VtkColorLookupTable::EXPONENTIAL);
			else
				lut->setInterpolationType(VtkColorLookupTable::NONE);
		}
		else // default
			lut->setInterpolationType(VtkColorLookupTable::NONE);

		QDomElement point = docElement.firstChildElement();
		double range[2] = { strtod((point.attribute("x")).toStdString().c_str(),0), strtod((point.attribute("x")).toStdString().c_str(),0) };

		while (!point.isNull())
		{
			if ((point.nodeName().compare("Point") == 0 )
				&& point.hasAttribute("x")
				&& point.hasAttribute("r")
				&& point.hasAttribute("g")
				&& point.hasAttribute("b"))
			{
				double value = strtod((point.attribute("x")).toStdString().c_str(),0);
				char r = static_cast<int>(255 * strtod((point.attribute("r")).toStdString().c_str(),0));
				char g = static_cast<int>(255 * strtod((point.attribute("g")).toStdString().c_str(),0));
				char b = static_cast<int>(255 * strtod((point.attribute("b")).toStdString().c_str(),0));
				char o = static_cast<int>(255 * (point.hasAttribute("o") ? strtod((point.attribute("o")).toStdString().c_str(),0) : 1));

				if (value<range[0]) range[0] = value;
				if (value>range[1]) range[1] = value;

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

#endif // XMLLUTREADER_H
