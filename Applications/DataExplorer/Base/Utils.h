/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <QListView>
#include <QStringList>
#include <QStringListModel>
#include <string>
#include <vector>

namespace Utils
{
std::vector<std::string> getSelectedObjects(QStringList const& list);

void moveSelectedItems(QListView* sourceView,
                       QStringListModel& sourceModel,
                       QStringListModel& targetModel);
}  // namespace Utils