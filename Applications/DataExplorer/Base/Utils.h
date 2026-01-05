// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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