/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-16
 * \brief  Definition of the StratScene class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QGraphicsScene>

#include "StationBorehole.h"
#include "Applications/DataHolderLib/Color.h"

class StratBar;
class QNonScalableGraphicsTextItem;

/**
 * \brief The scene for the visualisation of the stratigraphy of a borehole.
 */
class StratScene : public QGraphicsScene
{
public:
    /// Constructor
    StratScene(GeoLib::StationBorehole* station,
               std::map<std::string, DataHolderLib::Color>* stratColors = nullptr,
               QObject* parent = 0);
    ~StratScene();

    /// The margin between the boundary of the scene and the bounding box of all items within the scene
    static const int MARGIN = 50;

private:
    /// Adds text labels indicating the depth at the beginning and end of each soil layer
    void addDepthLabels(std::vector<GeoLib::Point*> profile, double offset);

    /// Add a non-scalable text item to the scene.
    QNonScalableGraphicsTextItem* addNonScalableText(const QString &text,
                                                     const QFont &font = QFont());

    /// Adds text labels indicating the name of each soil layer
    void addSoilNameLabels(std::vector<std::string> soilNames,
                           std::vector<GeoLib::Point*> profile,
                           double offset);

    /// Add a stratigraphy-bar to the scene.
    StratBar* addStratBar(GeoLib::StationBorehole* station,
                          std::map<std::string, DataHolderLib::Color>* stratColors = nullptr);
};
