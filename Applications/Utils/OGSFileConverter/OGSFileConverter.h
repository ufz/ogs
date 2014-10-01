/**
 * \file OGSFileConverter.h
 * 2012/04/04 KR Initial implementation
 */

#ifndef OGSFILECONVERTER_H
#define OGSFILECONVERTER_H

#include "ui_OGSFileConverter.h"
#include <QDialog>

class OGSFileConverter : public QDialog, private Ui_OGSFileConverter
{
	Q_OBJECT

public:
	OGSFileConverter(QWidget* parent = 0);
	~OGSFileConverter(void);

private:
	bool fileExists(const std::string &file_name) const;

private slots:
	void convertGML2GLI(const QStringList &input, const QString &output) const;
	void convertGLI2GML(const QStringList &input, const QString &output) const;
	void convertVTU2MSH(const QStringList &input, const QString &output) const;
	void convertMSH2VTU(const QStringList &input, const QString &output) const;

	void on_gml2gliButton_pressed() const;
	void on_gli2gmlButton_pressed() const;
	void on_vtu2mshButton_pressed() const;
	void on_msh2vtuButton_pressed() const;
	void on_closeDialogButton_pressed();
};

#endif //OGSFILECONVERTER_H
