/**
 * \file   SurfaceExtractionDialog.h
 * \author Karsten Rink
 * \date   2015-01-29
 * \brief  Definition of the SurfaceExtractionDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SURFACEEXTRACTIONDIALOG_H
#define SURFACEEXTRACTIONDIALOG_H

#include "ui_SurfaceExtraction.h"
#include <QDialog>

#include "MathLib/Vector3.h"

namespace MeshLib {
    class Mesh;
}

/**
 * \brief A dialog window for managing properties for writing meshes to files.
 */
class SurfaceExtractionDialog : public QDialog, private Ui_SurfaceExtraction
{
    Q_OBJECT

public:
    SurfaceExtractionDialog(QDialog* parent = 0);
    ~SurfaceExtractionDialog() {}

    int getTolerance() const { return _tolerance; }
    MathLib::Vector3 const& getNormal() const { return _dir; }

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject() { this->done(QDialog::Rejected); };

private:
    int _tolerance;
    MathLib::Vector3 _dir;
};

#endif //SURFACEEXTRACTIONDIALOG_H
