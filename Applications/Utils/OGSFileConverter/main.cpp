#include "OGSFileConverter.h"

#include <clocale>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include <QtGui/QApplication>

int main(int argc, char* argv[])
{
	ApplicationsLib::LogogSetup logog_setup;

	QApplication app(argc, argv);
	setlocale(LC_NUMERIC,"C");
	OGSFileConverter* fc = new OGSFileConverter();
	fc->setWindowTitle( fc->windowTitle() );
	fc->show();
	int returncode = app.exec();
	delete fc;

	return returncode;
}
