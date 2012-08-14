/**
 * \file OGSFileConverter.h
 * 2012/04/04 KR Initial implementation
 */

#ifndef OGSFILECONVERTER_H
#define OGSFILECONVERTER_H

#include "ui_OGSFileConverter.h"
#include <QDialog>

#include "FileFinder.h"

class OGSFileConverter : public QDialog, private Ui_OGSFileConverter
{
	Q_OBJECT

public:
	OGSFileConverter(QWidget* parent = 0);
	~OGSFileConverter(void);

private:
	FileFinder createFileFinder();

private slots:
	void convertGML2GLI(const QStringList &input, const QString &output);
	void convertGLI2GML(const QStringList &input, const QString &output);
	void convertCND2BC(const QStringList &input, const QString &output);
	void convertBC2CND(const QStringList &input, const QString &output);

	void on_gml2gliButton_pressed();
	void on_gli2gmlButton_pressed();
	void on_bc2cndButton_pressed();
	void on_cnd2bcButton_pressed();
	void on_closeDialogButton_pressed();

signals:

};

#endif //OGSFILECONVERTER_H
