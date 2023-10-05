/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <QStringList>
#include <string>
#include <vector>

namespace Utils
{
std::vector<std::string> getSelectedObjects(QStringList const& list);
}