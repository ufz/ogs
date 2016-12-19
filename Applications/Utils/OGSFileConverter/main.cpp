/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "OGSFileConverter.h"

#include <clocale>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include <QApplication>


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
