// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "QValueTooltipSlider.h"

#include <QCursor>
#include <QString>
#include <QToolTip>

QValueTooltipSlider::QValueTooltipSlider(QWidget* parent) : QSlider(parent)
{
    connect(this, SIGNAL(sliderMoved(int)), this, SLOT(setTooltipValue(int)));
}

void QValueTooltipSlider::setTooltipValue(int value)
{
    QToolTip::showText(QCursor::pos(), QString::number(value));
}
