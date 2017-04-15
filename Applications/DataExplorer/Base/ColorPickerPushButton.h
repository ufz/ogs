/**
 * \file
 * \author Lars Bilke
 * \date   2010-05-17
 * \brief  Definition of the ColorPickerPushButton class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    ColorPickerPushButton(QWidget* parent = nullptr);

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
