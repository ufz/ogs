/**
 * \file   OGSFileConverter.h
 * \author Karsten Rink
 * \date   2012-04-04
 * \brief  Definition of OGSFileConverter class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_OGSFileConverter.h"
#include <QDialog>

/**
 * A conversion tool for ogs5 and ogs6 files
 */
class OGSFileConverter : public QDialog, private Ui_OGSFileConverter
{
    Q_OBJECT

public:
    /// Constructor
    OGSFileConverter(QWidget* parent = nullptr);
    /// Destructor
    ~OGSFileConverter(void) override;

private:
    /// Checks if a given file already exists
    bool fileExists(const std::string &file_name) const;

private slots:
    /// Converts all files in the input list and writes the new files to the output directory with the same file name + updated extension.
    void convertGML2GLI(const QStringList &input, const QString &output) const;
    void convertGLI2GML(const QStringList &input, const QString &output) const;
    void convertVTU2MSH(const QStringList &input, const QString &output) const;
    void convertMSH2VTU(const QStringList &input, const QString &output) const;

    void on_gml2gliButton_pressed() const;
    void on_gli2gmlButton_pressed() const;
    void on_vtu2mshButton_pressed() const;
    void on_msh2vtuButton_pressed() const;
    void on_closeDialogButton_pressed();
};
