/**
 * \file QValueTooltipSlider.cpp
 * 23/03/2011 LB Initial implementation
 *
 * Implementation of QValueSlider class
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
