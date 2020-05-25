/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the DiagramList class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DiagramList.h"
#include "GetDateTime.h"

#include "BaseLib/Logging.h"

#include "DateTools.h"
#include "StringTools.h"
#include "SensorData.h"
#include <QFile>
#include <QTextStream>
#include <limits>

DiagramList::DiagramList() : xLabel_(""), yLabel_(""), xUnit_(""), yUnit_("") {}

DiagramList::~DiagramList() = default;

float DiagramList::calcMinXValue()
{
    auto min = std::min_element(
        coords_.begin(), coords_.end(),
        [](auto const& c0, auto const& c1) { return c0.first < c1.first; });
    if (min != coords_.end())
    {
        return min->first;
    }
    return std::numeric_limits<float>::max();
}

float DiagramList::calcMaxXValue()
{
    float max = std::numeric_limits<float>::lowest();
    std::size_t nCoords = coords_.size();
    for (std::size_t i = 0; i < nCoords; i++)
    {
        if (coords_[i].first > max)
        {
            max = coords_[i].first;
        }
    }
    return max;
}

float DiagramList::calcMinYValue()
{
    float min = std::numeric_limits<float>::max();
    std::size_t nCoords = coords_.size();
    for (std::size_t i = 0; i < nCoords; i++)
    {
        if (coords_[i].second < min)
        {
            min = coords_[i].second;
        }
    }
    return min;
}

float DiagramList::calcMaxYValue()
{
    float max = std::numeric_limits<float>::lowest();
    std::size_t nCoords = coords_.size();
    for (std::size_t i = 0; i < nCoords; i++)
    {
        if (coords_[i].second > max)
        {
            max = coords_[i].second;
        }
    }
    return max;
}

bool DiagramList::getPath(QPainterPath &path, float scaleX, float scaleY)
{
    QPointF p;
    if (getPoint(p, 0))
    {
        QPainterPath pp(QPointF(p.x() * scaleX, p.y() * scaleY));
        path = pp;

        std::size_t nCoords = coords_.size();
        for (std::size_t i = 1; i < nCoords; i++)
        {
            getPoint(p, i);
            path.lineTo(QPointF(p.x() * scaleX, p.y() * scaleY));
        }
        return true;
    }

    return false;
}

bool DiagramList::getPoint(QPointF &p, std::size_t i)
{
    if (i < coords_.size())
    {
        p.setX(coords_[i].first);
        p.setY(coords_[i].second);
        return true;
    }

    return false;
}

/*
 * Reads an external list into the coordinate arrays.
 * This method uses files containing the following format:
 *        xValue <tab> yValue
 * Both values may be int or double.
 */
/*
   int DiagramList::readList(char* path)
   {
    int date;
    double xVal, yVal;
    QString line;
    QStringList fields;

    QFile file(path);
    QTextStream in( &file );

    if (!file.open(QIODevice::ReadOnly))
    {
        return 0;
    }

    while (!in.atEnd()) {
        line = in.readLine();
        fields = line.split('\t');
        if (fields.size() >= 2) {
            xVal = fields.takeFirst().toDouble();
            yVal = fields.takeFirst().toDouble();
            xCoords.push_back(xVal);
            yCoords.push_back(yVal);
        }
        else return 0;
    }

    file.close();
    update();

    return 1;
   }*/

int DiagramList::readList(const QString &path, std::vector<DiagramList*> &lists)
{
    QFile file(path);
    QTextStream in( &file );

    if (!file.open(QIODevice::ReadOnly))
    {
        qDebug("Could not open file...");
        return 0;
    }

    QString line = in.readLine();
    QStringList fields = line.split('\t');
    int nLists(fields.size() - 1);

    if (fields.size() >= 2)
    {
        fields.takeFirst();
        for (int i = 0; i < nLists; i++)
        {
            auto* l = new DiagramList;
            l->setName(fields.takeFirst());
            //value = strtod(BaseLib::replaceStringreplaceString(",", ".", fields.takeFirst().toStdString()).c_str(),0);
            //l->setStartDate(startDate);
            //l->addNextPoint(0,value);
            lists.push_back(l);
        }

        bool first_loop(true);
        QDateTime startDate;
        QDateTime currentDate;
        unsigned line_count (1);

        while (!in.atEnd())
        {
            line = in.readLine();
            line_count++;
            fields = line.split('\t');
            if (fields.size() >= (nLists + 1))
            {
                QString const stringDate = fields.takeFirst();
                currentDate = getDateTime(stringDate);
                if (first_loop)
                {
                    startDate = currentDate;
                    for (int i = 0; i < nLists; i++)
                    {
                        lists[i]->setStartDate(startDate);
                    }
                    first_loop = false;
                }

                auto const numberOfSecs =
                    static_cast<float>(startDate.secsTo(currentDate));
                for (int i = 0; i < nLists; i++)
                {
                    float const value = static_cast<float>(
                        strtod(BaseLib::replaceString(
                                   ",", ".", fields.takeFirst().toStdString())
                                   .c_str(),
                               nullptr));
                    lists[i]->addNextPoint(numberOfSecs,value);
                }
            }
            else
            {
                WARN("DiagramList::readList(): Unexpected format in line {:d}.",
                     line_count);
                file.close();
                return 0;
            }
        }
    }
    else
    {
        qDebug("Unexpected file format...");
        file.close();
        return 0;
    }

    file.close();

    for (int i = 0; i < nLists; i++)
    {
        lists[i]->update();
    }

    return nLists;
}

int DiagramList::readList(const SensorData* data, std::vector<DiagramList*> &lists)
{
    std::vector<SensorDataType> const& time_series_names (data->getTimeSeriesNames());
    int nLists(time_series_names.size());

    std::vector<std::size_t> time_steps;
    if (data->getStepSize()>0)
    {
        const std::size_t start    = data->getStartTime();
        const std::size_t end      = data->getEndTime();
        const std::size_t stepsize = data->getStepSize();
        for (std::size_t i = start; i <= end; i += stepsize)
        {
            time_steps.push_back(i);
        }
    }
    else
    {
        time_steps = data->getTimeSteps();
    }

    bool is_date (false);

    if (!(BaseLib::int2date(time_steps[0])).empty())
    {
        is_date = true;
    }

    std::size_t nValues (time_steps.size());

    for (int i = 0; i < nLists; i++)
    {
        auto* l = new DiagramList;
        l->setName(QString::fromStdString(SensorData::convertSensorDataType2String(time_series_names[i])));
        l->setXLabel("Time");
        lists.push_back(l);

        const std::vector<float> *time_series = data->getTimeSeries(time_series_names[i]);

        if (is_date)
        {
            l->setXUnit("day");
            QDateTime const startDate(
                getDateTime(BaseLib::int2date(time_steps[0])));
            lists[i]->setStartDate(startDate);
            for (std::size_t j = 0; j < nValues; j++)
            {
                QDateTime const currentDate(
                    getDateTime(BaseLib::int2date(time_steps[j])));
                auto numberOfSecs =
                    static_cast<float>(startDate.secsTo(currentDate));
                lists[i]->addNextPoint(numberOfSecs, (*time_series)[j]);
            }
        }
        else
        {
            l->setXUnit("time step");
            for (std::size_t j = 0; j < nValues; j++)
            {
                lists[i]->addNextPoint(static_cast<float>(time_steps[j]),
                                       (*time_series)[j]);
            }
        }

        lists[i]->update();
    }

    return nLists;
}

void DiagramList::truncateToRange(QDateTime const& start, QDateTime const& end)
{
    auto start_secs = static_cast<float>(startDate_.secsTo(start));
    if (start_secs < 0)
    {
        start_secs = 0;
    }
    auto end_secs = static_cast<float>(startDate_.secsTo(end));
    if (end_secs < start_secs)
    {
        end_secs = coords_.back().first;
    }

    if (start_secs == 0 && end_secs == coords_.back().first)
    {
        return;
    }

    coords_.erase(
        std::remove_if(coords_.begin(), coords_.end(),
                       [&](std::pair<float, float> const& c) {
                           return (c.first < start_secs || c.first > end_secs);
                       }),
        coords_.end());
    startDate_ = start;
    for (auto& c : coords_)
    {
        c.first -= start_secs;
    }
    update();
}

void DiagramList::setList(
    std::vector<std::pair<QDateTime, float>> const& coords)
{
    if (coords.empty())
    {
        return;
    }

    startDate_ = coords[0].first;
    std::transform(coords.begin(), coords.end(), std::back_inserter(coords_),
                   [this](auto const& p) {
                       return std::make_pair(
                           static_cast<float>(startDate_.daysTo(p.first)),
                           p.second);
                   });

    update();
}

void DiagramList::setList(std::vector<std::pair<float, float>> const& coords)
{
    if (coords.empty())
    {
        return;
    }

    this->startDate_ = QDateTime();
    std::copy(coords.begin(), coords.end(), std::back_inserter(coords_));
    update();
}

std::size_t DiagramList::size() const
{
    return coords_.size();
}

void DiagramList::update()
{
    minX_ = calcMinXValue();
    maxX_ = calcMaxXValue();
    minY_ = calcMinYValue();
    maxY_ = calcMaxYValue();
}


