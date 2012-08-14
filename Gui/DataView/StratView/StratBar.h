/**
 * \file StratBar.h
 * 2010/03/16 - KR Initial implementation
 */

#ifndef STRATBAR_H
#define STRATBAR_H

#include "Station.h"
#include <QGraphicsItem>
#include <cmath>

/**
 * \brief A 2D bar visualisation of a borehole stratigraphy.
 *
 * A 2D bar visualisation of a borehole stratigraphy as a QGraphicsItem.
 */
class StratBar : public QGraphicsItem
{
public:
	/**
	 * \brief Constructor
	 * \param station The borehole whose stratigraphy will be visualised.
	 * \param stratColors A color map.
	 * \param parent The parent QGraphicsItem.
	 */
	StratBar(GEOLIB::StationBorehole* station,
	         std::map<std::string, GEOLIB::Color*>* stratColors = NULL,
	         QGraphicsItem* parent = 0);
	~StratBar();

	/// Returns the bounding rectangle of the bar.
	QRectF boundingRect() const;

	/// Paints the bar.
	void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);

private:
	/**
	 * \brief Calculates the height for a soil layer by "log(d+1)*100".
	 * \param h The original thickness of the soil layer.
	 */
	double logHeight(double h) const { return log(h + 1) * 100; }

	/// Calculates the total height of the bar by calculating and adding the log-height for all layers in the borehole
	double totalLogHeight() const;

	/// The default width of the bar
	static const int BARWIDTH = 50;

	GEOLIB::StationBorehole* _station;
	std::map<std::string, GEOLIB::Color*> _stratColors;
};

#endif //STRATBAR_H
