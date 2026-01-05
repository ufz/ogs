// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
