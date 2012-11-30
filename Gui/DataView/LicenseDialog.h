/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file LicenseDialog.h
 *
 * Created on 2012-11-30 by Karsten Rink
 */

#ifndef LICENSEDIALOG_H
#define LICENSEDIALOG_H

#include "ui_License.h"
#include <QDialog>

/**
 * \brief A dialog window displaying the OGS license information
 */
class LicenseDialog : public QDialog, private Ui_License
{
	Q_OBJECT

public:
	LicenseDialog(QDialog* parent = 0);
	~LicenseDialog() {};

private:
	void setText();

private slots:
	void on_okPushButton_pressed();

};

#endif //LICENSEDIALOG_H
