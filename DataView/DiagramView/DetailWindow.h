/**
 * \file DetailWindow.h
 * KR Initial implementation
 */

#ifndef DETAILWINDOW_H
#define DETAILWINDOW_H

#include <QtGui/QWidget>
#include "ui_DetailWindow.h"


/**
 * \brief Creates a window containing a diagram.
 */
class DetailWindow : public QWidget, private Ui_DetailWindow
{
	Q_OBJECT

public:
	/// Creates an empty diagram window.
	//DetailWindow(QWidget* parent = 0);
	/**
	 * Creates a window containing a diagram.
	 * \param filename ASCII file containing x and y values for the graph to be displayed.
	 * \param parent The parent QWidget.
	 */
	DetailWindow(QString filename, QWidget* parent = 0);

	/**
	 * Creates a window containing a diagram
	 * \param list A QDiagramList containing all the data points and necessary metainformation for a graph to be displayed
	 * \param parent The parent QWidget.
	 */
	DetailWindow(DiagramList* list, QWidget* parent = 0);
	~DetailWindow(void);

	/**
	 * Adds another plot to window. Axes are automatically resized, a random color is used.
	 */
	void addList(DiagramList* list);

	/**
	 * Adds another plot with a given colour to window. Axes are automatically resized.
	 */
	void addList(DiagramList* list, QColor c);

private:
	/// Automatically resize window based on the measurements of the included graphs.
	void resizeWindow();

	std::vector<DiagramList*> _list;

private slots:
	void on_closeButton_clicked();
};

#endif //DETAILWINDOW_H
