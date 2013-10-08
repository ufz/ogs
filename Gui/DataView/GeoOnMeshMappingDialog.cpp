/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-07
 * \brief  Implementation of the GeoOnMeshMappingDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoOnMeshMappingDialog.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "OGSError.h"


GeoOnMeshMappingDialog::GeoOnMeshMappingDialog(QDialog* parent)
	: _new_geo_name(""), QDialog(parent)
{
	setupUi(this);
}

GeoOnMeshMappingDialog::~GeoOnMeshMappingDialog()
{
}

void GeoOnMeshMappingDialog::accept()
{
	if (this->advancedMappingButton->isEnabled())
	{
		_new_geo_name = this->geoNameEdit->text().toStdString();
		if (_new_geo_name.empty())
		{
			OGSError::box("Please enter name for new geometry.");
			return;
		}
	}
	this->done(QDialog::Accepted);
}


