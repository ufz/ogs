// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QDialog>

namespace GeoLib {
struct GeoObject;
}

class QDialogButtonBox;
class QLabel;
class QComboBox;
class QVBoxLayout;

/**
 * \brief Small dialog for setting a name for an object.
 */
class SelectMeshDialog : public QDialog
{
    Q_OBJECT

public:
    /// Constructor
    SelectMeshDialog(const GeoLib::GeoObject* geo_object,
                     const std::list<std::string>& msh_names,
                     QDialog* parent = nullptr);
    ~SelectMeshDialog() override;

    QDialogButtonBox* _buttonBox; /// The buttons used in this dialog.

private:
    /// Constructs a dialog window
    void setupDialog(const std::list<std::string> &msh_names);

    QLabel* _txt_label;
    QComboBox* _msh_names;
    QVBoxLayout* _layout;
    const GeoLib::GeoObject* _geo_object;


private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

signals:
    //void requestNameChange(const std::string&, const GeoLib::GEOTYPE, std::size_t, std::string);
};
