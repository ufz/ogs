// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include <QApplication>
#include <clocale>

#include "BaseLib/Logging.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "OGSFileConverter.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "A conversion tool for ogs5 and ogs6 files.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> gmsh_path_arg(
        "g", "gmsh-path",
        "Input (.msh). The path to the input gmsh binary file", false, "",
        "INPUT_FILE");
    cmd.add(gmsh_path_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC, "C");
    auto* fc = new OGSFileConverter(gmsh_path_arg.getValue());
    fc->setWindowTitle(fc->windowTitle());
    fc->show();
    int returncode = QApplication::exec();
    delete fc;

    return returncode;
}
