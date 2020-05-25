/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the DiagramList class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QColor>
#include <QDateTime>
#include <QPainterPath>
#include <QPoint>
#include <vector>

class SensorData;

/**
 * \brief A List of data points and all the necessary meta-information to draw a graph.
 */
class DiagramList
{
public:
    /// Constructur containing an empty list.
    DiagramList();
    ~DiagramList();

    /// Returns the colour of the graph.
    QColor getColor() const {   return colour_; }

    /// Returns the height of the bounding box of all data points within the list.
    float height()    const { return maxY_ - minY_;  }

    /// Returns the minimum x-value.
    float minXValue() const { return minX_; }

    /// Returns the maximum x-value.
    float maxXValue() const { return maxX_; }

    /// Returns the minimum y-value.
    float minYValue() const { return minY_; }

    /// Returns the maximum y-value.
    float maxYValue() const { return maxY_; }

    /// Returns the start date of this list
    const QDateTime getStartDate() const { return startDate_; }

    /// Returns the name of the diagram.
    QString getName() const { return name_; }

    /**
     * Returns all the data points in form of a QPainterPath in scene coordinates.
     * The order of the points is the same as in the vector of coordinates.
     * \param path The path containing all data points.
     * \param scaleX Scaling factor in x-direction.
     * \param scaleY Scaling factor in y-direction.
     * \return true if everything is alright. false if the size of the coordinate arrays don't match.
     */
    bool getPath(QPainterPath &path, float scaleX, float scaleY);

    /**
     * Returns the position of one point in the vector of coordinates.
     * \param p The point-object that will be returned.
     * \param i Number of the point to be returned.
     * \return true if everything is alright. false if the point does not exist.
     */
    bool getPoint(QPointF &p, std::size_t i);

    /// Returns the label associated with the x-axis
    QString getXLabel() const { return xLabel_; }

    /// Returns the label associated with the y-axis
    QString getYLabel() const { return yLabel_; }

    /// Returns the unit associated with the x-axis
    QString getXUnit() const { return xUnit_; }

    /// Returns the unit associated with the y-axis
    QString getYUnit() const { return yUnit_; }

    /**
     * Reads information from a file. The reader assumes that values in the file are separated
     * by tabstops. Also, the first row should contain identifiers for the values and the first
     * column should contain timestamps. Currently accepted timestamps are of the following
     * formats:
     *   "dd.mm.yyyy"
     *   "dd.mm.yyyy.hh.mm.ss" (this is the timestamp format used for the UFZ database)
     * For each column after the timestamps a new diagram list is created.
     */
    static int readList(const QString &path, std::vector<DiagramList*> &list);

    static int readList(const SensorData* data, std::vector<DiagramList*> &list);

    /// Sets the colour of the graph.
    void setColor(QColor c) { colour_ = c; }

    /// Sets the name of the graph to be displayed in the caption.
    void setName(QString name) { name_ = name; }

    /// Adds a point at (x,y) to the list
    void addNextPoint(float x, float y) { coords_.emplace_back(x, y); }

    /// cut list entries not within the given range
    void truncateToRange(QDateTime const& start, QDateTime const& end);

    /// Sets the start date (i.e. the min-value of the x-axis).
    void setStartDate(QDateTime date) { startDate_ = date; }

    /// Specifies the meaning of the x Axis.
    void setXLabel(QString label) { xLabel_ = label; }

    /// Specifies the meaning of the y Axis.
    void setYLabel(QString label) { yLabel_ = label; }

    /**
     * Sets the list of x/y-coordinates.
     * \param coords List of coordinates.
     */
    void setList(std::vector<std::pair<float, float>> const& coords);

    /**
     * Sets the list of x/y-coordinates.
     * Note: This function converts QDateTime values to float values of the number of
     * days from the first date (which is set as day 0)
     * \param coords List of coordinates.
     */
    void setList(std::vector<std::pair<QDateTime, float>> const& coords);

    /// Specifies the unit of the x Axis.
    void setXUnit(QString unit) { xUnit_ = unit; }

    /// Specifies the unit of the y Axis.
    void setYUnit(QString unit) { yUnit_ = unit; }

    /// Returns the number of data points.
    std::size_t size() const;

    /// Returns the width of the bounding box of all data points within the list.
    double width() const { return maxX_ - minX_; }

private:
    /// Returns the minimum x-value of all the data points.
    float calcMinXValue();

    /// Returns the maximum x-value of all the data points.
    float calcMaxXValue();

    /// Returns the minimum y-value of all the data points.
    float calcMinYValue();

    /// Returns the maximum y-value of all the data points.
    float calcMaxYValue();

    /**
     * Reads an ASCII file containing the coordinates in the following format:
     *        date (tab) value
     * where 'date' is given as 'dd.mm.yyyy'.
     * (Changes to that format are easily implemented using QTimeDate)
     * \return Returns 1 if everything is alright. Returns 0 and displays an error message if there was an error.
     */
    int readLine(std::ifstream inFile, QDateTime &cDate, float &cValue);

    /// Updates the bounds of the data points contained in the list.
    void update();

    float maxX_{0};
    float maxY_{0};
    float minX_{0};
    float minY_{0};
    std::vector< std::pair<float, float> > coords_;
    QString name_;
    QString xLabel_;
    QString yLabel_;
    QString xUnit_;
    QString yUnit_;
    QColor colour_;
    QDateTime startDate_;
};
