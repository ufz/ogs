/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file QValueTooltipSlider.h
 *
 * Created on 2011-03-23 by Lars Bilke
 */

#ifndef QVALUETOOLTIPSLIDER_H
#define QVALUETOOLTIPSLIDER_H

#include <QSlider>

/**
 * \brief A QSlider which shows its value as a tooltip when moved.
 */
class QValueTooltipSlider : public QSlider
{
	Q_OBJECT

public:
	QValueTooltipSlider(QWidget* parent = 0);

public slots:
	void setTooltipValue(int value);

protected:
};

#endif // QVALUETOOLTIPSLIDER_H
