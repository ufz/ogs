// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QGraphicsScene>

#include "GeoLib/StationBorehole.h"
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
    explicit StratScene(
        GeoLib::StationBorehole* station,
        std::map<std::string, DataHolderLib::Color>* stratColors = nullptr,
        QObject* parent = nullptr);
    ~StratScene() override;

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
