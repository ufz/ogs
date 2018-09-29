/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-16
 * \brief  Definition of the StratView class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "StratScene.h"
#include <QGraphicsView>
#include <QWidget>


namespace GeoLib
{
class StationBorehole;
}

/**
 * \brief A view in which to display the stratigraphy of a borehole.
 */
class StratView : public QGraphicsView
{
    Q_OBJECT
public:
    /**
     * Creates an empty view.
     */
    StratView(QWidget* parent = nullptr) : _scene(nullptr) { Q_UNUSED(parent); }
    ~StratView() override;

    /// Sets the Borehole whose data should be visualised.
    void setStation(GeoLib::StationBorehole* station,
                    std::map<std::string, DataHolderLib::Color>* stratColors = nullptr);

    /// Returns the height of the bounding rectangle of all objects within the scene.
    int getHeight() { return static_cast<int>((_scene->itemsBoundingRect()).height()); }

    /// Returns the width of the bounding rectangle of all objects within the scene.
    int getWidth() { return static_cast<int>((_scene->itemsBoundingRect()).width()); }

    void saveAsImage(QString fileName);

protected:
    /// Resizes the scene.
    void resizeEvent(QResizeEvent* event) override;

private:
    /// Initialises the view.
    void initialize();

    /// The minimum size of the window.
    QSize minimumSizeHint() const override
    {
        return QSize(3 * _scene->MARGIN, 2 * _scene->MARGIN);
    }

    /// The default size of the window.
    QSize sizeHint() const override
    {
        return QSize(6 * _scene->MARGIN, 4 * _scene->MARGIN);
    }

    /// Updates the view automatically when a Borehole is added or when the window containing the view changes its state.
    void update();

    StratScene* _scene;
};
