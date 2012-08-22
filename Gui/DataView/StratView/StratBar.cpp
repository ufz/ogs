/**
 * \file StratBar.cpp
 * 2010/03/16 - KR Initial implementation
 */

#include "StratBar.h"
#include <QPainter>

StratBar::StratBar(GeoLib::StationBorehole* station,
                   std::map<std::string, GeoLib::Color*>* stratColors,
                   QGraphicsItem* parent) :
	QGraphicsItem(parent), _station(station)
{
	if (stratColors)
		_stratColors = *stratColors;
}

StratBar::~StratBar()
{
//	std::map<std::string, GeoLib::Color*>::iterator it;
//	for (it = _stratColors.begin(); it != _stratColors.end(); it++) {
//		delete it->second;
//	}
}

QRectF StratBar::boundingRect() const
{
	return QRectF(0, 0, BARWIDTH + 10, totalLogHeight());
}

void StratBar::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
	Q_UNUSED (option)
	Q_UNUSED (widget)

	double top = 0, height = 0;

	QPen pen(Qt::black, 1, Qt::SolidLine, Qt::SquareCap, Qt::RoundJoin);
	pen.setCosmetic(true);
	painter->setPen(pen);
	//painter->drawRect(_bar);

	//pen.setWidth(1);
	std::vector<GeoLib::Point*> profile = _station->getProfile();
	std::vector<std::string> soilNames = _station->getSoilNames();
	size_t nLayers = profile.size();

	painter->drawLine(0, 0, BARWIDTH + 5, 0);

	for (size_t i = 1; i < nLayers; i++)
	{
		top += height;
		height = logHeight(((*(profile[i - 1]))[2] - (*(profile[i]))[2]));
		QRectF layer(0, top, BARWIDTH, height);
		const GeoLib::Color* c (GeoLib::getColor(soilNames[i], _stratColors));
		QBrush brush(QColor((int)(*c)[0],
		                    (int)(*c)[1],
		                    (int)(*c)[2],
		                    127), Qt::SolidPattern);
		painter->setBrush(brush);

		painter->drawRect(layer);
		painter->drawLine(0, (int)layer.bottom(), BARWIDTH + 5, (int)layer.bottom());
		//painter->drawText(BARWIDTH+10, layer.bottom(), QString::number((*(profile[i]))[2]));
	}
}

double StratBar::totalLogHeight() const
{
	double height = 0;
	std::vector<GeoLib::Point*> profile = _station->getProfile();

	for (size_t i = 1; i < profile.size(); i++)
		height += ( log((*(profile[i - 1]))[2] - (*(profile[i]))[2] + 1) * 100 );

	return height;
}
