/**
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "OGSFileConverter.h"

#include <clocale>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "BaseLib/BuildInfo.h"

#include <QApplication>


int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "A conversion tool for ogs5 and ogs6 files.\n\n"
        "OpenGeoSys-6 software, version " +
            BaseLib::BuildInfo::git_describe +
            ".\n"
            "Copyright (c) 2012-2019, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> gmsh_path_arg("g", "gmsh-path",
                                               "the path to the gmsh binary",
                                               false, "", "path as string");
    cmd.add(gmsh_path_arg);
    cmd.parse( argc, argv );
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");
    auto* fc = new OGSFileConverter(gmsh_path_arg.getValue());
    fc->setWindowTitle( fc->windowTitle() );
    fc->show();
    int returncode = app.exec();
    delete fc;

    return returncode;
}
