/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-17
 * \brief  Definition of the ColorTableModel class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QAbstractTableModel>
#include <QColor>

#include "Applications/DataHolderLib/Color.h"

/**
 * The PolylinesModel is a Qt model which represents Polylines.
 */
class ColorTableModel : public QAbstractTableModel
{
    Q_OBJECT

public:
    ColorTableModel(
        const std::map<std::string, DataHolderLib::Color*>& colorLookupTable,
        QObject* parent = nullptr);
    ~ColorTableModel() override;

    int columnCount(const QModelIndex& parent = QModelIndex()) const override;

    QVariant data(const QModelIndex& index, int role) const override;

    int rowCount(const QModelIndex& parent = QModelIndex()) const override
    {
        Q_UNUSED (parent);
        return _listOfPairs.size();
    }

    QVariant headerData(int section, Qt::Orientation orientation,
                        int role /*= Qt::DisplayRole*/) const override;

private:
    bool buildTable( const std::map<std::string, DataHolderLib::Color*> &colorLookupTable );

    QList< QPair<QString, QColor> > _listOfPairs;
};
