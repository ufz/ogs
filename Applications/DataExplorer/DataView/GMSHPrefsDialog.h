// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ui_GMSHPrefs.h"
#include <QDialog>

class QStringListModel;

namespace GeoLib
{
class GEOObjects;
}

/**
 * \brief A dialog window for setting preferences for GMSH
 */
class GMSHPrefsDialog : public QDialog, private Ui_GMSHPrefs
{
    Q_OBJECT

public:
    explicit GMSHPrefsDialog(GeoLib::GEOObjects const& geoObjects,
                             QDialog* parent = nullptr);
    ~GMSHPrefsDialog() override;

private:
    std::vector<std::string> getSelectedObjects(QStringList list);

    QStringListModel* _allGeo;
    QStringListModel* _selGeo;

private slots:
    void on_selectGeoButton_pressed();
    void on_deselectGeoButton_pressed();
    void on_radioAdaptive_toggled(bool isTrue);

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

signals:
    void requestMeshing(std::vector<std::string> &, unsigned, double, double, double, bool);
};
