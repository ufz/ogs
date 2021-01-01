/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "OGSFileConverter.h"

#include <clocale>
#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"
#include <QApplication>


int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "A conversion tool for ogs5 and ogs6 files.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
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
    int returncode = QApplication::exec();
    delete fc;

    return returncode;
}
