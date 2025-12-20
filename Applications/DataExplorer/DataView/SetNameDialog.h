// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QDialog>
#include "GeoLib/GeoType.h"

class QDialogButtonBox;
class QLabel;
class QLineEdit;
class QVBoxLayout;

/**
 * \brief Small dialog for setting a name for an object.
 */
class SetNameDialog : public QDialog
{
    Q_OBJECT

public:
    /// Constructor
    SetNameDialog(const std::string& geo_object_type, std::size_t id,
                  const std::string& old_name, QDialog* parent = nullptr);
    ~SetNameDialog() override;

    std::string getNewName();

private:
    /// Constructs a dialog window
    void setupDialog(const QString &label, const std::string &old_name);

    QLabel* _txt_label;
    QLineEdit* _new_name;
    QVBoxLayout* _layout;
    QDialogButtonBox* _buttonBox;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;
};
