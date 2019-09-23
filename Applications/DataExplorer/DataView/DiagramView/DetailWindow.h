/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the DetailWindow class.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_DetailWindow.h"
#include <QWidget>


/**
 * \brief Creates a window containing a diagram.
 */
class DetailWindow : public QWidget, private Ui_DetailWindow
{
    Q_OBJECT

public:
    /// Creates an empty diagram window.
    explicit DetailWindow(QWidget* parent = nullptr);
    /**
     * Creates a window containing a diagram.
     * \param filename ASCII file containing x and y values for the graph to be displayed.
     * \param parent The parent QWidget.
     */
    explicit DetailWindow(QString filename, QWidget* parent = nullptr);

    /**
     * Creates a window containing a diagram
     * \param list A QDiagramList containing all the data points and necessary metainformation for a graph to be displayed
     * \param parent The parent QWidget.
     */
    explicit DetailWindow(DiagramList* list, QWidget* parent = nullptr);

    explicit DetailWindow(std::vector<std::size_t> data,
                          QWidget* parent = nullptr);

    ~DetailWindow() override;

    /**
     * Adds another plot to window. Axes are automatically resized, a random color is used.
     */
    void addList(DiagramList* list);

    /**
     * Adds another plot with a given colour to window. Axes are automatically resized.
     */
    void addList(DiagramList* list, QColor c);

private:
    /// Automatically resize window based on the measurements of the included graphs.
    void resizeWindow();

private slots:
    void on_addDataButton_clicked();
    void on_closeButton_clicked();
};
