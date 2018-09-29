/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-17
 * \brief  Implementation of the ColorTableView class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "ColorTableView.h"
#include <QHeaderView>
#include <QPainter>

ColorTableView::ColorTableView( QWidget* parent /*= 0*/ ) : QTableView(parent)
{
    this->verticalHeader()->hide();
    this->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->resizeColumnsToContents();
    this->resizeRowsToContents();
}

void ColorTableViewDelegate::paint(QPainter* painter,
                                   const QStyleOptionViewItem &option,
                                   const QModelIndex &index) const
{
    QColor val;
    if (index.column() == 1)
    {
        if (index.data().canConvert(QMetaType::QColor))
        {
            val = index.data().value<QColor>();
            QBrush brush(val);
            painter->fillRect(option.rect, brush);
        }
    }
    else
        QItemDelegate::paint(painter, option, index);
}

QSize ColorTableViewDelegate::sizeHint( const QStyleOptionViewItem &option,
                                        const QModelIndex &index ) const
{
    QSize s = QItemDelegate::sizeHint(option, index);
    if( s.isValid() )
        s.setHeight((int)(0.5 * s.height()));
    return s;
}
