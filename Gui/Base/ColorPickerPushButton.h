/**
 * \file ColorPickerPushButton.h
 * 17/5/2010 LB Initial implementation
 *
 */

#ifndef COLORPICKERPUSHBUTTON_H
#define COLORPICKERPUSHBUTTON_H

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
	ColorPickerPushButton(QWidget* parent = 0);

public slots:
	/// Calls the QColorDialog
	void mouseReleaseEvent(QMouseEvent* e);

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

#endif // COLORPICKERPUSHBUTTON_H
