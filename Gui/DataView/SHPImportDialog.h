/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file SHPImportDialog.h
 *
 * Created on 2010-01-25 by Karsten Rink
 */

#ifndef SHPIMPORTDIALOG_H
#define SHPIMPORTDIALOG_H

#include <QDialog>

class SHPInterface;
class GEOModels;

class QDialogButtonBox;
class QFileInfo;
class QGridLayout;
class QLabel;
class QLineEdit;
class QRadioButton;
class QVBoxLayout;

/**
 * \brief Dialog for selecting which information should be loaded from a shape file.
 */
class SHPImportDialog : public QDialog
{
	Q_OBJECT

public:
	/// Constructor
	SHPImportDialog(std::string filename, GEOModels* geoModels, QDialog* parent = 0);
	~SHPImportDialog();

	QDialogButtonBox* _buttonBox; /// The buttons used in this dialog.

private:
	/// Constructs a dialog window based on the information found in the selected shape file
	void setupDialog();

	QGridLayout* _layout;
	QLabel* _shpContentLabel;
	QLabel* _nameLabel;
	QLineEdit* _listName;
	QRadioButton* _choice1, * _choice2;
	std::string _filename;
	short _fileType;
	SHPInterface* _shpInterface;

private slots:
	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void shpLoaded(QString);
};

#endif //SHPIMPORTDIALOG_H
