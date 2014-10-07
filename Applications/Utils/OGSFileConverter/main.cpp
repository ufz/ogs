#include "OGSFileConverter.h"

#include "logog/include/logog.hpp"
#include "LogogSimpleFormatter.h"

#include <QtGui/QApplication>

int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;
	BaseLib::LogogSimpleFormatter* formatter = new BaseLib::LogogSimpleFormatter;
	logogCout->SetFormatter(*formatter);

	QApplication app(argc, argv);
	setlocale(LC_NUMERIC,"C");
	OGSFileConverter* fc = new OGSFileConverter();
	fc->setWindowTitle( fc->windowTitle() );
	fc->show();
	int returncode = app.exec();
	delete fc;

	delete formatter;
	delete logogCout;
	LOGOG_SHUTDOWN();

	return returncode;
}
