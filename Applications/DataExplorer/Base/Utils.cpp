/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Base/Utils.h"

namespace Utils
{
std::vector<std::string> getSelectedObjects(QStringList const& list)
{
    std::vector<std::string> indexList;
    std::transform(list.begin(), list.end(), std::back_inserter(indexList),
                   [](auto const& index) { return index.toStdString(); });
    return indexList;
}

void moveSelectedItems(QListView* sourceView,
                       QStringListModel& sourceModel,
                       QStringListModel& targetModel)
{
    QModelIndexList selected = sourceView->selectionModel()->selectedIndexes();
    QStringList targetList = targetModel.stringList();

    for (const auto& index : selected)
    {
        targetList.append(index.data().toString());
        sourceModel.removeRow(index.row());
    }

    targetModel.setStringList(targetList);
}
}  // namespace Utils