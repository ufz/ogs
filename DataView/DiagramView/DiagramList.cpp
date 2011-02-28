/**
 * \file DiagramList.cpp
 * KR Initial implementation
 */

#include <limits>
#include <QFile>
#include <QTextStream>
#include "DiagramList.h"
#include "StringTools.h"


DiagramList::DiagramList() : _maxX(0), _maxY(0), _minX(0), _minY(0), _xLabel(""), _yLabel(""), _xUnit(""), _yUnit(""), _startDate()
{
}

DiagramList::~DiagramList()
{
}

float DiagramList::calcMinXValue()
{
	float min = std::numeric_limits<float>::max();
	size_t nCoords = _coords.size();
	for (size_t i=0; i<nCoords; i++)
	{
		if (_coords[i].first<min) min=_coords[i].first;
	}
	return min;
}

float DiagramList::calcMaxXValue()
{
	float max = std::numeric_limits<float>::min();
	size_t nCoords = _coords.size();
	for (size_t i=0; i<nCoords; i++)
	{
		if (_coords[i].first>max) max=_coords[i].first;
	}
	return max;
}

float DiagramList::calcMinYValue()
{
	float min = std::numeric_limits<float>::max();
	size_t nCoords = _coords.size();
	for (size_t i=0; i<nCoords; i++)
	{
		if (_coords[i].second<min) min=_coords[i].second;
	}
	return min;
}

float DiagramList::calcMaxYValue()
{
	float max = std::numeric_limits<float>::min();
	size_t nCoords = _coords.size();
	for (size_t i=0; i<nCoords; i++)
	{
		if (_coords[i].second>max) max=_coords[i].second;
	}
	return max;
}

bool DiagramList::getPath(QPainterPath &path, float scaleX, float scaleY)
{
	QPointF p;
	if (getPoint(p,0))
	{
		QPainterPath pp(QPointF(p.x()*scaleX, p.y()*scaleY));
		path = pp;

		size_t nCoords = _coords.size();
		for (size_t i=1; i<nCoords; i++)
		{
			getPoint(p,i);
			path.lineTo(QPointF(p.x()*scaleX, p.y()*scaleY));
		}
		return true;
	}
	else
		return false;
}

bool DiagramList::getPoint(QPointF &p, size_t i)
{
	if (i<_coords.size())
	{
		p.setX(_coords[i].first);
		p.setY(_coords[i].second);
		return true;
	}
	else return false;
}


/*
 * Reads an external list into the coordinate arrays.
 * This method uses files containing the following format:
 *		xValue <tab> yValue
 * Both values may be int or double.
 */
/*
int DiagramList::readList(char* path)
{
	int date;
	double xVal, yVal;
	QString line;
	QStringList fields;

	QFile file(path);
	QTextStream in( &file );

	if (!file.open(QIODevice::ReadOnly))
    {
		return 0;
	}

	while (!in.atEnd()) {
		line = in.readLine();
		fields = line.split('\t');
		if (fields.size() >= 2) {
			xVal = fields.takeFirst().toDouble();
			yVal = fields.takeFirst().toDouble();
			xCoords.push_back(xVal);
			yCoords.push_back(yVal);
		}
		else return 0;
	}

    file.close();
	update();

    return 1;
}*/

int DiagramList::readList(const QString &path, std::vector<DiagramList*> &list)
{
	QFile file(path);
	QTextStream in( &file );

	if (!file.open(QIODevice::ReadOnly))
    {
		qDebug("Could not open file...");
		return 0;
	}

	QString line = in.readLine();
	QStringList fields = line.split('\t');
	int nLists(fields.size()-1);

	
	if (fields.size() >= 2) 
	{
		int numberOfDays(0);
		QString stringDate(fields.takeFirst());
		double value(0);
		QDateTime startDate(QDateTime::fromString(stringDate, "dd.MM.yyyy"));

		for (int i=0; i<nLists; i++)
		{
			DiagramList* l = new DiagramList;
			value = strtod(replaceString(",", ".", fields.takeFirst().toStdString()).c_str(),0);
			l->setStartDate(startDate);
			l->addNextPoint(0,value);
			list.push_back(l);
		}

		QDateTime currentDate;

		while (!in.atEnd()) {
			line = in.readLine();
			fields = line.split('\t');
			if (fields.size() >= (nLists+1)) {
				stringDate = fields.takeFirst();
				currentDate = QDateTime::fromString(stringDate, "dd.MM.yyyy");
				numberOfDays = startDate.daysTo(currentDate);

				for (int i=0; i<nLists; i++)
				{
					value = strtod(replaceString(",", ".", fields.takeFirst().toStdString()).c_str(),0);
					list[i]->addNextPoint(numberOfDays,value);
				}
			}
			else 
			{
				qDebug("Unexpected file format...");
				file.close();
				return 0;
			}
		}
	}
	else
	{
		qDebug("Unexpected file format...");
		file.close();
		return 0;
	}

	file.close();
	
	for (int i=0; i<nLists; i++)
		list[i]->update();
	
    return nLists;
}

void DiagramList::setList(std::vector< std::pair<QDateTime, float> > coords)
{
	int numberOfDays;
	QDateTime startDate;

	startDate = coords[0].first;
	_coords.push_back(std::pair<float, float>(0, coords[0].second));

	size_t nCoords = coords.size();
	for (size_t i=1; i<nCoords; i++)
	{
		numberOfDays = startDate.daysTo(coords[i].first);
		_coords.push_back(std::pair<float, float>(numberOfDays, coords[i].second));
	}

	update();
}

void DiagramList::setList(std::vector< std::pair<float, float> > coords)
{
	size_t nCoords = coords.size();
	for (size_t i=0; i<nCoords; i++)
		_coords.push_back(coords[i]);

	update();
}

size_t DiagramList::size()
{
	if (!(_coords.empty()))
		return _coords.size();
	else return 0;
}

void DiagramList::update()
{
	_minX = calcMinXValue();
	_maxX = calcMaxXValue();
	_minY = calcMinYValue();
	_maxY = calcMaxYValue();
}
