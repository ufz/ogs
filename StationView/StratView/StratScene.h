/**
 * \file StratScene.h
 * 2010/03/16 - KR Initial implementation
 */

#ifndef STRATSCENE_H
#define STRATSCENE_H

#include <QGraphicsScene>
#include "Station.h"

class StratBar;
class QNonScalableGraphicsTextItem;

/**
 * \brief The scene for the visualisation of the stratigraphy of a borehole.
 */
class StratScene : public QGraphicsScene
{
public:
	/// Constructor
	StratScene(GEOLIB::StationBorehole* station, std::map<std::string, GEOLIB::Color*> *stratColors = NULL, QObject* parent = 0);
	~StratScene();

	/// The margin between the boundary of the scene and the bounding box of all items within the scene
	static const int MARGIN=50;

private:
	/// Adds text labels indicating the depth at the beginning and end of each soil layer
	void addDepthLabels(std::vector<GEOLIB::Point*> profile, double offset);

	/// Add a non-scalable text item to the scene.
	QNonScalableGraphicsTextItem* addNonScalableText(const QString &text, const QFont &font = QFont());

	/// Adds text labels indicating the name of each soil layer
	void addSoilNameLabels(std::vector<std::string> soilNames, std::vector<GEOLIB::Point*> profile, double offset);

	/// Add a stratigraphy-bar to the scene.
	StratBar* addStratBar(GEOLIB::StationBorehole* station, std::map<std::string, GEOLIB::Color*> *stratColors = NULL);

};

#endif //STRATSCENE_H
