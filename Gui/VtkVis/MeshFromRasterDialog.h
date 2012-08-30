/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MeshFromRasterDialog.h
 *
 * Created on 2011-11-17 by Karsten Rink
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
