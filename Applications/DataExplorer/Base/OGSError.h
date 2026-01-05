// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

class QString;

/**
 * \brief Manages error messages via message boxes
 */
class OGSError
{
public:
    /**
     * Displays an error in a QMessageBox
     * \param e The error message.
     */
    static void box(const QString &e);

    /**
     * Displays an error in a QMessageBox
     * \param e The error message.
     * \param t The title of the message box
     * \sa QMessageBox
     */
    static void box(const QString &e, const QString &t);

    /**
     * Displays a question in a QMessageBox (offering Ok | Cancel options)
     * Default value is 'Cancel' so that no bad things happen if the user
     * presses enter without reading the text (e.g. when overwriting files)
     * \param e The error message.
     * \param t The title of the message box
     * \return 'true' if 'Ok' has been pressed, 'false' otherwise
     * \sa QMessageBox
     */
    static bool question(const QString &e, const QString &t);


protected:
    OGSError();
    ~OGSError();
};
