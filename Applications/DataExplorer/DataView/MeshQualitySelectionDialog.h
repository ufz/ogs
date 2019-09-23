/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-16
 * \brief  Definition of the MshQualitySelectionDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_MeshQualitySelection.h"

#include <QDialog>
#include "MeshLib/MeshEnums.h"

class VtkMeshSource;

/**
 * \brief A dialog for selecting a mesh quality metric
 */
class MeshQualitySelectionDialog : public QDialog, private Ui_MeshQualitySelection
{
    Q_OBJECT

public:
    explicit MeshQualitySelectionDialog(QDialog* parent = nullptr);
    ~MeshQualitySelectionDialog() override;

    /// Returns selected metric
    MeshLib::MeshQualityType getSelectedMetric() const { return _metric; }

    /// Returns true if a histogram needs to be calculated
    bool getHistogram() const { return this->histogramCheckBox->isChecked(); }

    /// Returns selected path for histogram (or empty string if no histogram is required)
    std::string getHistogramPath() const { return _histogram_path; }

private:
    MeshLib::MeshQualityType _metric{MeshLib::MeshQualityType::EDGERATIO};
    std::string _histogram_path;

private slots:
    void on_histogramCheckBox_toggled(bool is_checked) const;
    void on_histogramPathButton_pressed();

    void accept() override;
    void reject() override;
};
