// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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