/**
 * \file CondFromRasterDialog.h
 * 2012/01/04 KR Initial implementation
 */

#ifndef CONDFROMRASTERDIALOG_H
#define CONDFROMRASTERDIALOG_H

#include "ui_CondFromRaster.h"
#include <QDialog>

#include "ProjectData.h"


/**
 * \brief A dialog window for creating DIRECT boundary conditions from raster files
 */
class CondFromRasterDialog : public QDialog, private Ui_CondFromRaster
{
	Q_OBJECT

public:
	CondFromRasterDialog(const ProjectData &project, QDialog* parent = 0);
	~CondFromRasterDialog(void);

private:
	const ProjectData _project;

private slots:
	void on_selectButton_pressed();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();
	
signals:

};

#endif //CONDFROMRASTERDIALOG_H
