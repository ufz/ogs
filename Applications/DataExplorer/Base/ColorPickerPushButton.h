// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// ** INCLUDES **
#include <QPushButton>

class QColor;
class QMouseEvent;

/**
 * ColorPickerPushButton calls a QColorDialog on clicking and then sends a
 * colorPicked(QColor) signal. It also saves the last color and sets its
 * background color accordingly.
 */
class ColorPickerPushButton : public QPushButton
{
    Q_OBJECT

public:
    explicit ColorPickerPushButton(QWidget* parent = nullptr);

public slots:
    /// Calls the QColorDialog
    void mouseReleaseEvent(QMouseEvent* e) override;

    /// Sets the color.
    void setColor(QColor color);
    void setColor(double* color);

private:
    QString colorToCss(QColor color);
    QString colorToString(QColor color);

    QColor _color;

signals:
    /// Is emitted when a color was picked from the dialog.
    void colorPicked(QColor);
};
