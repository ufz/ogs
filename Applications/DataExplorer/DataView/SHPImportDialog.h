/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-25
 * \brief  Definition of the SHPImportDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QDialog>

namespace FileIO {
    class SHPInterface;
}

namespace GeoLib {
    class GEOObjects;
}

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
    SHPImportDialog(std::string filename, GeoLib::GEOObjects& geo_objects,
                    std::string const& gmsh_path,
                    QDialog* parent = nullptr);
    ~SHPImportDialog() override;

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
    FileIO::SHPInterface* _shpInterface;
    std::string const _gmsh_path;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

signals:
    void shpLoaded(QString);
};
