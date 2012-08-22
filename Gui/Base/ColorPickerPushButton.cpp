/**
 * \file ColorPickerPushButton.cpp
 * 17/5/2010 LB Initial implementation
 *
 * Implementation of ColorPickerPushButton
 */

// ** INCLUDES **
#include "ColorPickerPushButton.h"

#include <QColorDialog>

ColorPickerPushButton::ColorPickerPushButton( QWidget* parent /*= 0*/ )
	: QPushButton(parent)
{
	setAutoFillBackground(true);
	_color = QColor("white");
}

void ColorPickerPushButton::mouseReleaseEvent(QMouseEvent* e)
{
	Q_UNUSED(e);
	QColor newColor = QColorDialog::getColor(_color, NULL, "Choose a color");
	if (!newColor.isValid())
		return;

	setColor(newColor);

	emit colorPicked(_color);
}

QString ColorPickerPushButton::colorToCss( QColor color )
{
	QString colorStr = "rgb";
	colorStr.append(colorToString(color));

	return colorStr;
}

QString ColorPickerPushButton::colorToString( QColor color )
{
	QString colorStr = "(";
	colorStr.append(QString::number(color.red()));
	colorStr.append(", ");
	colorStr.append(QString::number(color.green()));
	colorStr.append(", ");
	colorStr.append(QString::number(color.blue()));
	colorStr.append(")");

	return colorStr;
}

void ColorPickerPushButton::setColor( QColor color )
{
	_color = color;

	// Compute text color
	QColor hsv = _color.toHsv();
	QString textColorStr;
	if (hsv.valueF() < 0.5f)
		textColorStr = "color: rgb(255, 255, 255);";
	else
		textColorStr = "color: rgb(0, 0, 0);";

	QString stylesheetStr = "background-color: ";
	stylesheetStr.append(colorToCss(_color));
	stylesheetStr.append(";");
	stylesheetStr.append(textColorStr);
	this->setStyleSheet(stylesheetStr);

	this->setText(colorToString(_color));
}

void ColorPickerPushButton::setColor( double* color )
{
	return setColor(QColor::fromRgbF(color[0], color[1], color[2]));
}
