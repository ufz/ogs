/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-07
 * \brief  Definition of the GeoOnMeshMappingDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOTOMESHMAPPINGDIALOG_H
#define GEOTOMESHMAPPINGDIALOG_H

#include "ui_GeoOnMeshMapping.h"
#include <QDialog>

#include "GeoMapper.h"

/**
 * \brief A dialog window for creating DIRECT boundary conditions from raster files
 */
class GeoOnMeshMappingDialog : public QDialog, private Ui_GeoOnMeshMapping
{
	Q_OBJECT

public:
	GeoOnMeshMappingDialog(QDialog* parent = 0);
	~GeoOnMeshMappingDialog(void);

	std::string getNewGeoName() const { return _new_geo_name; };

private:
	std::string _new_geo_name;

private slots:
	void on_advancedMappingButton_toggled(bool isSelected) { this->geoNameEdit->setEnabled(isSelected); };

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject() { this->done(QDialog::Rejected); };

};

#endif //GEOTOMESHMAPPINGDIALOG_H
