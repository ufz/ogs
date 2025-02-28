/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-16
 * \brief  Definition of the StratBar class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cmath>

#include <QGraphicsItem>

#include "GeoLib/StationBorehole.h"
#include "Applications/DataHolderLib/Color.h"

/**
 * \brief A 2D bar visualisation of a borehole stratigraphy.
 *
 * A 2D bar visualisation of a borehole stratigraphy as a QGraphicsItem.
 */
class StratBar : public QGraphicsItem
{
public:
    /**
     * \brief Constructor
     * \param station The borehole whose stratigraphy will be visualised.
     * \param stratColors A color map.
     * \param parent The parent QGraphicsItem.
     */
    explicit StratBar(
        GeoLib::StationBorehole* station,
        std::map<std::string, DataHolderLib::Color>* stratColors = nullptr,
        QGraphicsItem* parent = nullptr);
    ~StratBar() override = default;

    /// Returns the bounding rectangle of the bar.
    QRectF boundingRect() const override;

    /// Paints the bar.
    void paint(QPainter* painter,
               const QStyleOptionGraphicsItem* option,
               QWidget* widget) override;

private:
    /**
     * \brief Calculates the height for a soil layer by "log(d+1)*100".
     * \param h The original thickness of the soil layer.
     */
    double logHeight(double h) const { return log(h + 1) * 100; }

    /// Calculates the total height of the bar by calculating and adding the log-height for all layers in the borehole
    double totalLogHeight() const;

    /// The default width of the bar
    static const int BARWIDTH = 50;

    GeoLib::StationBorehole* _station;
    std::map<std::string, DataHolderLib::Color> _stratColors;
};
