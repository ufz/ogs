/**
 * \file VisPrefsDialog.h
 * 14/06/2010 KR Initial implementation
 */

#ifndef VISPREFSDIALOG_H
#define VISPREFSDIALOG_H

#include <QtGui/QMainWindow>
#include "ui_VisPrefs.h"
#include "Point.h"

class VtkVisPipeline;

/**
 * \brief A dialog window for settung up a visualisation preferences
 */
class VisPrefsDialog : public QDialog, private Ui_VisPrefsDialog
{
	Q_OBJECT

public:
	VisPrefsDialog(VtkVisPipeline* pipeline, QDialog* parent = 0);
	~VisPrefsDialog(void);




private slots:
	/// Sets the background colour.
	void on_bgColorButton_colorPicked(QColor color);

	/// Adds a light above the scene.
	void on_lightAboveBox_clicked();

	/// Adds a light below the scene.
	void on_lightBelowBox_clicked();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

private:
	VtkVisPipeline* _vtkVisPipeline;
	GEOLIB::Point _above;
	GEOLIB::Point _below;

};

#endif //VISPREFSDIALOG_H
