/**
 * \file MeshFromRasterDialog.h
 * 2011/11/17 KR Initial implementation
 */

#ifndef MSHFROMRASTERDIALOG_H
#define MSHFROMRASTERDIALOG_H

#include "ui_MeshFromRaster.h"
#include "VtkMeshConverter.h"

#include <QtGui/QDialog>


/**
 * \brief A dialog for specifying the parameters to construct a mesh based on a raster
 */
class MeshFromRasterDialog : public QDialog, private Ui_MeshFromRaster
{
	Q_OBJECT

public:
	/// Constructor
	MeshFromRasterDialog(QDialog* parent = 0);

	~MeshFromRasterDialog(void);

private slots:
	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void setMeshParameters(QString, MshElemType::type, UseIntensityAs::type);

};

#endif //MSHFROMRASTERDIALOG_H
