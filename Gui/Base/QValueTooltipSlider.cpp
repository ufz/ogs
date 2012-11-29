/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file QValueTooltipSlider.cpp
 *
 * Created on 2011-03-23 by Lars Bilke
 */

// ** INCLUDES **
#include "QValueTooltipSlider.h"

#include <QCursor>
#include <QString>
#include <QToolTip>

#include <iostream>

QValueTooltipSlider::QValueTooltipSlider(QWidget* parent)
	: QSlider(parent)
{
	connect(this, SIGNAL(sliderMoved(int)), this, SLOT(setTooltipValue(int)));
}

void QValueTooltipSlider::setTooltipValue(int value)
{
	QToolTip::showText(QCursor::pos(), QString::number(value));
}
