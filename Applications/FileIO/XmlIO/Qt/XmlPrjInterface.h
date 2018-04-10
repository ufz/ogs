/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <QString>
#include <QtXml/QDomNode>

#include "BaseLib/IO/XmlIO/Qt/XMLQtInterface.h"
#include "BaseLib/IO/XmlIO/XMLInterface.h"

#include "Applications/DataHolderLib/BoundaryCondition.h"
#include "Applications/DataHolderLib/Project.h"
#include "Applications/DataHolderLib/SourceTerm.h"

namespace FileIO
{
/**
 * Data Explorer XML interface for project files
 */
class XmlPrjInterface : public BaseLib::IO::XMLInterface,
                        public BaseLib::IO::XMLQtInterface
{
public:
    XmlPrjInterface(DataHolderLib::Project& project);

    ~XmlPrjInterface() override = default;

    /// Reads an xml-file containing a project.
    int readFile(const QString& fileName) override;

    /// Reads an xml-file containing a project.
    bool readFile(std::string const& fname) override
    {
        return readFile(QString(fname.c_str())) != 0;
    }

    /// Writes a project to the specified file
    int writeToFile(const std::string& filename);

protected:
    bool write() override;

private:
    /// Tests if a given parameter exists within the file
    QDomNode findParam(QDomNode const& node, QString const& param_name) const;

    /// Manages reading all kinds of conditions
    void readConditions(QDomNode const& node, QDomNode const& param_root);

    /// Reading all boundary conditions
    void readBoundaryConditions(QDomNode const& bc_root,
                                QDomNode const& param_root,
                                DataHolderLib::ProcessVariable const& pvar);

    /// Reading all source terms
    void readSourceTerms(QDomNode const& st_root, QDomNode const& param_root,
                         DataHolderLib::ProcessVariable const& pvar);

    /// Writes information on process variables
    void writeProcessVariables(QDomDocument& doc, QDomElement& root) const;

    /// Compiles a vector of all existing primary variables for writing purposes
    std::vector<DataHolderLib::ProcessVariable> getPrimaryVariableVec() const;

    /// Writes one specific condition
    template <typename T>
    void writeCondition(QDomDocument& doc, QDomElement& tag,
                        DataHolderLib::FemCondition const& cond) const;

    /// Writes a list of boundary conditions
    void writeBoundaryConditions(QDomDocument& doc, QDomElement& bc_list_tag,
                                 std::string const& name) const;

    /// Writes a list of source terms
    void writeSourceTerms(QDomDocument& doc, QDomElement& bc_list_tag,
                          std::string const& name) const;

    /// Parsing one specific condition
    template <typename T>
    T* parseCondition(QDomNode const& st_root, QDomNode const& param_root,
                      DataHolderLib::ProcessVariable const& pvar) const;

    std::string _filename;
    DataHolderLib::Project& _project;
};

}  // end namespace FileIO
