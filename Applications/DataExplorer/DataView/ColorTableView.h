// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QItemDelegate>
#include <QTableView>

/**
 *    A QTableView to display colour lookup tables.
 */
class ColorTableView : public QTableView
{
    Q_OBJECT

public:
    /// Constructor
    explicit ColorTableView(QWidget* parent = nullptr);
};

/**
 *    A delegate class to manage properties of ColorTableView.
 */
class ColorTableViewDelegate : public QItemDelegate
{
    Q_OBJECT

public:
    /// Constructor
    explicit ColorTableViewDelegate(QWidget* parent = nullptr)
        : QItemDelegate(parent)
    {
    }
    /// Overwrites the paint-method to set user-defined properties instead of the default properties.
    void paint(QPainter* painter, const QStyleOptionViewItem& option,
               const QModelIndex& index) const override;

    QSize sizeHint(const QStyleOptionViewItem& option,
                   const QModelIndex& index) const override;
};
