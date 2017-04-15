/**
 * \file
 * \author Lars Bilke
 * \date   2010-08-19
 * \brief  Implementation of the CheckboxDelegate class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "CheckboxDelegate.h"
#include <QApplication>
#include <QCheckBox>
#include <QEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QStyleOptionButton>

CheckboxDelegate::CheckboxDelegate(QObject* parent)
    : QItemDelegate(parent)
{
}

void CheckboxDelegate::paint(QPainter* painter, const QStyleOptionViewItem& option,
                             const QModelIndex& index) const
{
    if(index.isValid())
    {
        bool checked = index.model()->data(index, Qt::DisplayRole).toBool();

        QStyleOptionButton styleOptionButton;
        styleOptionButton.state |= QStyle::State_Enabled;
        if (checked)
            styleOptionButton.state |= QStyle::State_On;
        else
            styleOptionButton.state |= QStyle::State_Off;

        styleOptionButton.rect = this->checkboxRect(option);

        QApplication::style()->drawControl(QStyle::CE_CheckBox,
                                           &styleOptionButton, painter);
    }
    else
        QItemDelegate::paint(painter, option, index);
}

bool CheckboxDelegate::editorEvent(QEvent* event, QAbstractItemModel* model,
                                   const QStyleOptionViewItem &option, const QModelIndex &index)
{
    Q_UNUSED(option);

    if ((event->type() == QEvent::MouseButtonRelease) ||
        (event->type() == QEvent::MouseButtonDblClick))
    {
        auto* mouse_event = static_cast<QMouseEvent*>(event);
        if (mouse_event->button() != Qt::LeftButton ||
            !checkboxRect(option).contains(mouse_event->pos()))
            return false;
        if (event->type() == QEvent::MouseButtonDblClick)
            return true;
    }
    else if (event->type() == QEvent::KeyPress)
    {
        if (static_cast<QKeyEvent*>(event)->key() != Qt::Key_Space &&
            static_cast<QKeyEvent*>(event)->key() != Qt::Key_Select)
            return false;
    }
    else
        return false;

    bool checked = index.model()->data(index, Qt::DisplayRole).toBool();
    return model->setData(index, !checked, Qt::EditRole);
}

QSize CheckboxDelegate::sizeHint(const QStyleOptionViewItem & option,
                                 const QModelIndex & index) const
{
    Q_UNUSED(index);

    QRect rect = checkboxRect(option);
    return QSize(rect.width(), rect.height());
}

QRect CheckboxDelegate::checkboxRect(const QStyleOptionViewItem& viewItemStyleOptions) const
{
    QStyleOptionButton styleOptionButton;
    QRect rect = QApplication::style()->subElementRect(
            QStyle::SE_CheckBoxIndicator, &styleOptionButton);
    QPoint point(viewItemStyleOptions.rect.x() +
                 viewItemStyleOptions.rect.width() / 2 - rect.width() / 2,
                 viewItemStyleOptions.rect.y() + viewItemStyleOptions.rect.height() / 2 -
                 rect.height() / 2);
    return QRect(point, rect.size());
}
