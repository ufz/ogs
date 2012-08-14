#include "OGSFileConverter.h"
#include <QtGui/QApplication>

int main(int argc, char* argv[])
{
	QApplication app(argc, argv);
	setlocale(LC_NUMERIC,"C");
	OGSFileConverter* fc = new OGSFileConverter();
	fc->setWindowTitle( fc->windowTitle() );
	fc->show();
	int returncode = app.exec();
	delete fc;

	return returncode;
}
