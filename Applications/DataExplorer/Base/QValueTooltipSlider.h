// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QSlider>

/**
 * \brief A QSlider which shows its value as a tooltip when moved.
 */
class QValueTooltipSlider : public QSlider
{
    Q_OBJECT

public:
    explicit QValueTooltipSlider(QWidget* parent = nullptr);

public slots:
    void setTooltipValue(int value);

protected:
};
