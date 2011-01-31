/**
 * \file StratView.h
 * 2010/03/16 - KR Initial implementation
 */

#ifndef STRATVIEW_H
#define STRATVIEW_H

#include <QtGui/QWidget>
#include <QGraphicsView>
#include "StratScene.h"

namespace GEOLIB
{
	class StationBorehole;
}

/**
 * \brief A view in which to display the stratigraphy of a borehole.
 */
class StratView : public QGraphicsView
{
	Q_OBJECT
public:
	/**
	 * Creates an empty view.
	 */
	StratView(QWidget* parent = 0) : _scene(NULL) {Q_UNUSED(parent);}
	~StratView();

	/// Sets the Borehole whose data should be visualised.
	void setStation(GEOLIB::StationBorehole* station, std::map<std::string, GEOLIB::Color*> *stratColors = NULL);

	/// Returns the height of the bounding rectangle of all objects within the scene.
	int getHeight() { return static_cast<int>((_scene->itemsBoundingRect()).height()); }

	/// Returns the width of the bounding rectangle of all objects within the scene.
	int getWidth() { return static_cast<int>((_scene->itemsBoundingRect()).width()); }

	void saveAsImage(QString fileName);

protected:
	/// Resizes the scene.
	void resizeEvent(QResizeEvent* event);

private:
    /// Initialises the view.
	void initialize();

	/// The minimum size of the window.
	QSize minimumSizeHint() const { return QSize(3*_scene->MARGIN,2*_scene->MARGIN); }

	/// The default size of the window.
	QSize sizeHint() const { return QSize(6*_scene->MARGIN, 4*_scene->MARGIN); }

	/// Updates the view automatically when a Borehole is added or when the window containing the view changes its state.
	void update();

	StratScene* _scene;
};

#endif //STRATVIEW_H
