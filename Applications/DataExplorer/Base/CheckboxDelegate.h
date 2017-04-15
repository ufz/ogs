/**
 * \file
 * \author Lars Bilke
 * \date   2010-08-19
 * \brief  Definition of the CheckboxDelegate class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QItemDelegate>

class QWidget;
class QRect;

/**
 * \brief CheckboxDelegate modifies a model view to display boolean values as checkboxes.
 *
 * Important: the column on which this delegate is set (QAbstractItemView::setItemDelegateForColumn())
 * must not have the flags Qt::ItemIsEditable or Qt::ItemIsUserCheckable set in the model.
 **/
class CheckboxDelegate : public QItemDelegate
{
    Q_OBJECT

public:
    /// \brief Constructor
    CheckboxDelegate(QObject* parent = nullptr);

    /// \brief Paints a checkbox. This overrides the default painting of a combo box.
    void paint(QPainter* painter, const QStyleOptionViewItem& option,
               const QModelIndex& index) const;

    /// \brief Handles the click events and sets the model data.
    bool editorEvent(QEvent* event, QAbstractItemModel* model,
                     const QStyleOptionViewItem &option, const QModelIndex &index);

    QSize sizeHint (const QStyleOptionViewItem & option, const QModelIndex & index) const;

private:
    QRect checkboxRect(const QStyleOptionViewItem& viewItemStyleOptions) const;
};
