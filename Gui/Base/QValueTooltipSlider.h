/**
 * \file QValueTooltipSlider.h
 * 23/03/2011 LB Initial implementation
 */

#ifndef QVALUETOOLTIPSLIDER_H
#define QVALUETOOLTIPSLIDER_H

#include <QSlider>

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
