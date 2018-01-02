/**
 * \file   XmlNumInterface.cpp
 * \author Karsten Rink
 * \date   2014-08-05
 * \brief  Implementation of the XmlNumInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "XmlNumInterface.h"

#include <QFile>
#include <QTextCodec>
#include <QtXml/QDomDocument>

#include <logog/include/logog.hpp>

#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileFinder.h"


namespace FileIO
{
XmlNumInterface::XmlNumInterface() :
XMLInterface(), XMLQtInterface(BaseLib::FileFinder({BaseLib::BuildInfo::app_xml_schema_path}).getPath("OpenGeoSysNUM.xsd"))
{
}

int XmlNumInterface::readFile(QString const& fileName)
{
    if(XMLQtInterface::readFile(fileName) == 0)
        return 0;

    QDomDocument doc("OGS-NUM-DOM");
    doc.setContent(_fileData);
    QDomElement const docElement = doc.documentElement(); //OGSNonlinearSolverSetup
    if (docElement.nodeName().compare("OGSNonlinearSolverSetup"))
    {
        ERR("XmlNumInterface::readFile() - Unexpected XML root.");
        return 0;
    }

    QDomElement num_node = docElement.firstChildElement();

    while (!num_node.isNull())
    {
        if (num_node.nodeName().compare("Type") == 0)
        {
            std::string const solver_type = num_node.toElement().text().toStdString();
            INFO("Non-linear solver type: %s.", solver_type.c_str());
        }
        else if (num_node.nodeName().compare("LinearSolver") == 0)
            readLinearSolverConfiguration(num_node);
        else if (num_node.nodeName().compare("IterationScheme") == 0)
            readIterationScheme(num_node);
        else if (num_node.nodeName().compare("Convergence") == 0)
            readConvergenceCriteria(num_node);

        num_node = num_node.nextSiblingElement();
    }

    return 1;
}

void XmlNumInterface::readLinearSolverConfiguration(QDomElement const& lin_root)
{
    std::string const library = lin_root.attribute("Library").toStdString();
    std::string lin_solver_type, precond_type("no");

    QDomElement linear_solver_node = lin_root.firstChildElement();
    while (!linear_solver_node.isNull())
    {
        if (linear_solver_node.nodeName().compare("Type") == 0)
            lin_solver_type = linear_solver_node.toElement().text().toStdString();
        if (linear_solver_node.nodeName().compare("Preconditioner") == 0)
            precond_type = linear_solver_node.toElement().text().toStdString();
        linear_solver_node = linear_solver_node.nextSiblingElement();
    }
    INFO("Using %s-library with solver %s and %s preconditioner.", library.c_str(), lin_solver_type.c_str(), precond_type.c_str());
}


void XmlNumInterface::readIterationScheme(QDomElement const& iteration_root)
{
    QDomElement iteration_node = iteration_root.firstChildElement();
    int max_iterations(0), fixed_step_size(0);

    while (!iteration_node.isNull())
    {
        if (iteration_node.nodeName().compare("MaxIterations") == 0)
            max_iterations = iteration_node.toElement().text().toInt();
        if (iteration_node.nodeName().compare("FixedStepSize") == 0)
            fixed_step_size = iteration_node.toElement().text().toDouble();

        iteration_node = iteration_node.nextSiblingElement();
    }
    INFO("Doing a maximum of %d iterations at fixed step size of %f", max_iterations, fixed_step_size);
}

void XmlNumInterface::readConvergenceCriteria(QDomElement const& convergence_root)
{
    QDomElement conv_node = convergence_root.firstChildElement();
    double error_threshold(0);
    std::string error_method("");

    while (!conv_node.isNull())
    {
        if (conv_node.nodeName().compare("Method") == 0)
            error_method = conv_node.toElement().text().toStdString();
        if (conv_node.nodeName().compare("ErrorThreshold") == 0)
            error_threshold = conv_node.toElement().text().toDouble();

        conv_node = conv_node.nextSiblingElement();
    }
    INFO("Convergence reached when error below %f using %s.", error_threshold, error_method.c_str());
}

bool XmlNumInterface::write()
{
    INFO("Not yet implemented...");
    return false;
}

} //namespace
