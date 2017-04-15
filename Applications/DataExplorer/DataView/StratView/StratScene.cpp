/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-16
 * \brief  Implementation of the StratScene class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <limits>
#include <math.h>

#include <QGraphicsTextItem>

#include "DateTools.h"
#include "QNonScalableGraphicsTextItem.h"
#include "StratBar.h"
#include "StratScene.h"

StratScene::StratScene(GeoLib::StationBorehole* station,
                       std::map<std::string, DataHolderLib::Color>* stratColors,
                       QObject* parent) : QGraphicsScene(parent)
{
    QRectF textBounds;
    int stratBarOffset = 250;

    QFont font( "Arial", 15, QFont::DemiBold, false);

    QNonScalableGraphicsTextItem* boreholeTag = addNonScalableText("Borehole", font);
    QNonScalableGraphicsTextItem* boreholeName = addNonScalableText(
            "\"" + QString::fromStdString(station->getName()) + "\"",
            font);
    textBounds = boreholeTag->boundingRect();
    boreholeTag->setPos((textBounds.width() / 2.0), 80);
    textBounds = boreholeName->boundingRect();
    boreholeName->setPos((textBounds.width() / 2.0), 200);

    QNonScalableGraphicsTextItem* totalDepth =
            addNonScalableText("Depth: " + QString::number(station->getDepth()) + " m");
    textBounds = totalDepth->boundingRect();
    totalDepth->setPos((textBounds.width() / 2.0), 350);
/*
    QNonScalableGraphicsTextItem* dateText = addNonScalableText("Date: " + QString::fromStdString(date2string(station->getDate())));
    textBounds = dateText->boundingRect();
    dateText->setPos(this->MARGIN + (textBounds.width()/2.0), 350);
 */
    QNonScalableGraphicsTextItem* dot = addNonScalableText(" ");
    dot->setPos(0, 0);

    StratBar* stratBar = addStratBar(station, stratColors);
    stratBar->setPos(stratBarOffset, MARGIN);
    QRectF stratBarBounds = stratBar->boundingRect();

    addDepthLabels(station->getProfile(), stratBarOffset + stratBarBounds.width());

    if (station->getSoilNames().size() > 0)
        addSoilNameLabels(station->getSoilNames(), station->getProfile(), stratBarOffset +
                          (stratBarBounds.width() / 2));
}

StratScene::~StratScene() = default;

void StratScene::addDepthLabels(std::vector<GeoLib::Point*> profile, double offset)
{
    QRectF textBounds;
    double vertPos = MARGIN;
    std::vector<QNonScalableGraphicsTextItem*> depthText;
    depthText.push_back(addNonScalableText(QString::number((*(profile[0]))[2])));
    textBounds = depthText[0]->boundingRect();
    depthText[0]->setPos(offset + textBounds.width() / 2, vertPos);

    for (std::size_t i = 1; i < profile.size(); i++)
    {
        depthText.push_back(addNonScalableText(QString::number((*(profile[i]))[2])));
        vertPos += log((*(profile[i - 1]))[2] - (*(profile[i]))[2] + 1) * 100;
        textBounds = depthText[i]->boundingRect();
        depthText[i]->setPos(offset + textBounds.width() / 2, vertPos);
    }
}

QNonScalableGraphicsTextItem* StratScene::addNonScalableText(const QString &text, const QFont &font)
{
    auto* item = new QNonScalableGraphicsTextItem(text);
    item->setFont(font);
    addItem(item);
    return item;
}

void StratScene::addSoilNameLabels(std::vector<std::string> soilNames,
                                   std::vector<GeoLib::Point*> profile,
                                   double offset)
{
    //QRectF textBounds;
    double vertPos = MARGIN, halfHeight = 0;
    std::vector<QNonScalableGraphicsTextItem*> soilText;
    soilText.push_back(addNonScalableText(QString::fromStdString(soilNames[0])));
    //textBounds = soilText[0]->boundingRect();
    soilText[0]->setPos(offset /* - textBounds.width() */, vertPos);

    for (std::size_t i = 1; i < soilNames.size(); i++)
    {
        soilText.push_back(addNonScalableText(QString::fromStdString(soilNames[i])));
        halfHeight = log((*(profile[i - 1]))[2] - (*(profile[i]))[2] + 1) * 100 / 2;
        //textBounds = soilText[i]->boundingRect();
        soilText[i]->setPos(offset /* - textBounds.width() */, vertPos + halfHeight);
        vertPos += ( 2 * halfHeight );
    }
}

StratBar* StratScene::addStratBar(GeoLib::StationBorehole* station,
                                  std::map<std::string, DataHolderLib::Color>* stratColors)
{
    auto* b = new StratBar(station, stratColors);
    addItem(b);
    return b;
}
