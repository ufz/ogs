/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-17
 * \brief  Definition of the MeshFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MSHFROMRASTERDIALOG_H
#define MSHFROMRASTERDIALOG_H

#include "ui_MeshFromRaster.h"

#include <string>
#include <QDialog>

namespace MeshLib {
    enum class MeshElemType;
    enum class UseIntensityAs;
}


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

    std::string getMeshName() const { return _mesh_name; }
    std::string getArrayName() const { return _array_name; }
    MeshLib::MeshElemType getElementSelection() const { return _element_selection; }
    MeshLib::UseIntensityAs getIntensitySelection() const { return _intensity_selection; }

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject();

private:
    std::string _mesh_name;
    std::string _array_name;
    MeshLib::MeshElemType _element_selection;
    MeshLib::UseIntensityAs _intensity_selection;

};

#endif //MSHFROMRASTERDIALOG_H
