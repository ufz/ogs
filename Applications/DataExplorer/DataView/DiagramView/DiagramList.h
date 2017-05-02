/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the DiagramList class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    QColor getColor() const {   return _colour; }

    /// Returns the height of the bounding box of all data points within the list.
    float height()    const { return _maxY - _minY;  }

    /// Returns the minimum x-value.
    float minXValue() const { return _minX; }

    /// Returns the maximum x-value.
    float maxXValue() const { return _maxX; }

    /// Returns the minimum y-value.
    float minYValue() const { return _minY; }

    /// Returns the maximum y-value.
    float maxYValue() const { return _maxY; }

    /// Returns the start date of this list
    const QDateTime getStartDate() const { return _startDate; }

    /// Returns the name of the diagram.
    QString getName() const { return _name; }

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
    QString getXLabel() const { return _xLabel; }

    /// Returns the label associated with the y-axis
    QString getYLabel() const { return _yLabel; }

    /// Returns the unit associated with the x-axis
    QString getXUnit() const { return _xUnit; }

    /// Returns the unit associated with the y-axis
    QString getYUnit() const { return _yUnit; }

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
    void setColor(QColor c) { _colour = c; }

    /// Sets the name of the graph to be displayed in the caption.
    void setName(QString name) { _name = name; }

    /// Adds a point at (x,y) to the list
    void addNextPoint(float x, float y) { _coords.emplace_back(x, y); }
    /// Sets the start date (i.e. the min-value of the x-axis).
    void setStartDate(QDateTime date) { _startDate = date; }

    /// Specifies the meaning of the x Axis.
    void setXLabel(QString label) { _xLabel = label; }

    /// Specifies the meaning of the y Axis.
    void setYLabel(QString label) { _yLabel = label; }

    /**
     * Sets the list of x/y-coordinates.
     * \param coords List of coordinates.
     */
    void setList(std::vector< std::pair<float, float> > coords);

    /**
     * Sets the list of x/y-coordinates.
     * Note: This function converts QDateTime values to float values of the number of
     * days from the first date (which is set as day 0)
     * \param coords List of coordinates.
     */
    void setList(std::vector< std::pair<QDateTime, float> > coords);

    /// Specifies the unit of the x Axis.
    void setXUnit(QString unit) { _xUnit = unit; }

    /// Specifies the unit of the y Axis.
    void setYUnit(QString unit) { _yUnit = unit; }

    /// Returns the number of data points.
    std::size_t size();

    /// Returns the width of the bounding box of all data points within the list.
    double width() const { return _maxX - _minX; }

private:
    /// Returns the minimum x-value of all the data points.
    float calcMinXValue();

    /// Returns the maximum x-value of all the data points.
    float calcMaxXValue();

    /// Returns the minimum y-value of all the data points.
    float calcMinYValue();

    /// Returns the maximum y-value of all the data points.
    float calcMaxYValue();

    static QDateTime getDateTime(QString s);

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

    float _maxX;
    float _maxY;
    float _minX;
    float _minY;
    std::vector< std::pair<float, float> > _coords;
    QString _name;
    QString _xLabel;
    QString _yLabel;
    QString _xUnit;
    QString _yUnit;
    QColor _colour;
    QDateTime _startDate;
};
