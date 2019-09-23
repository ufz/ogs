/**
 * \file
 * \author Karsten Rink
 * \date   2014-08-05
 * \brief  Definition of the XmlNumInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/IO/XmlIO/XMLInterface.h"
#include "BaseLib/IO/XmlIO/Qt/XMLQtInterface.h"

namespace FileIO
{

class XmlNumInterface : public BaseLib::IO::XMLInterface, public BaseLib::IO::XMLQtInterface
{
public:
    XmlNumInterface();

    ~XmlNumInterface() override = default;

    int readFile(QString const& fileName) override;

    bool readFile(std::string const& fname) override
    {
        return readFile(QString(fname.c_str())) != 0;
    }

protected:
    void readLinearSolverConfiguration(QDomElement const& lin_root);
    void readIterationScheme(QDomElement const& iteration_root);
    void readConvergenceCriteria(QDomElement const& convergence_root);
    bool write() override;
};

}  // namespace FileIO
