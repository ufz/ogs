// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"

#include <vector>

class vtkThreshold;

class VtkColorLookupTable;

/// @brief This filter selects/thresholds elements based on the selected array.
class VtkCompositeElementSelectionFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeElementSelectionFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeElementSelectionFilter() override = default;

    void init() override;

    // Sets the range for the quality measure (default is [0,1] but this may vary for area- and volume-metrics).
    void setRange(double min_val, double max_val) { _range = std::make_pair(min_val, max_val); }

    void setSelectionArray(const std::string &selection_name, const std::vector<double> &selection = std::vector<double>());

    void SetUserVectorProperty(QString name, QList<QVariant> values) override;

private:
    /// Returns a colour lookup table optimised for quality measures
    VtkColorLookupTable* GetLookupTable();

    std::pair<double, double> _range;
    std::string _selection_name;
    std::vector<double> _selection;
};
