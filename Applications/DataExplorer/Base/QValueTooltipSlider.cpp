/**
 * \file
 * \author Lars Bilke
 * \date   2011-03-23
 * \brief  Implementation of the QValueTooltipSlider class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "QValueTooltipSlider.h"

#include <QCursor>
#include <QString>
#include <QToolTip>

QValueTooltipSlider::QValueTooltipSlider(QWidget* parent)
    : QSlider(parent)
{
    connect(this, SIGNAL(sliderMoved(int)), this, SLOT(setTooltipValue(int)));
}

void QValueTooltipSlider::setTooltipValue(int value)
{
    QToolTip::showText(QCursor::pos(), QString::number(value));
}
