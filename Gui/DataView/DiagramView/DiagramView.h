/**
 * \file DiagramView.h
 * KR Initial implementation
 */

#ifndef DIAGRAMVIEW_H
#define DIAGRAMVIEW_H

#include "DiagramScene.h"
#include <QGraphicsView>
#include <QtGui/QWidget>

/**
 * \brief A view in which to display a DiagramScene.
 *
 * A view in which to display a DiagramScene. It supports resizing of the window and loading of data into the diagram.
 */
class DiagramView : public QGraphicsView
{
	Q_OBJECT
public:
	/**
	 * Creates an empty view.
	 */
	DiagramView(QWidget* parent = 0);
	/**
	 * Creates a view already containing a graph
	 * \param list Contains a list of data points and metainformation to be displayed by the scene.
	 * \param parent The parent QWidget.
	 */
	DiagramView(DiagramList* list, QWidget* parent = 0);
	~DiagramView();

	/// Adds a new graph to the scene.
	void addGraph(DiagramList* list);
	/// Returns the height of the bounding rectangle of all objects within the scene.
	int getHeight();
	/// Returns the width of the bounding rectangle of all objects within the scene.
	int getWidth();

protected:
	/// Resizes the scene.
	void resizeEvent(QResizeEvent* event);

private:
	void initialize();
	void keepItemAspectRatio();
	QSize minimumSizeHint() const;
	QSize sizeHint() const;
	void update();

	DiagramScene* _scene;
};

#endif //DIAGRAMVIEW_H
