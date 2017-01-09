/**
 * \file
 * \author Karsten Rink
 * \date   2012-11-30
 * \brief  Implementation of the LicenseDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LicenseDialog.h"
#include <QTextStream>

LicenseDialog::LicenseDialog(QDialog* parent) : QDialog(parent)
{
    setupUi(this);
    this->licenseTextBrowser->setOpenExternalLinks(true);
    this->setText();
}

void LicenseDialog::setText()
{
    QString text;
    QTextStream stream(&text);
    stream << "<p>Copyright (c) 2012-2017, OpenGeoSys Community "
           << "(<a href='http://www.opengeosys.org'>http://www.opengeosys.org</a>)<br />"
           << "All rights reserved.</p>"
           << "<p>Redistribution and use in source and binary forms, with or without"
           << "modification, are permitted provided that the following conditions are met:"
           << "<ol>"
           << "<li>Redistributions of source code must retain the above copyright"
              << "notice, this list of conditions and the following disclaimer.</li>"
              << "<li>Redistributions in binary form must reproduce the above copyright"
              << "notice, this list of conditions and the following disclaimer in the"
              << "documentation and/or other materials provided with the distribution.</li>"
              << "<li>All advertising materials mentioning features or use of this software"
              << "must display the following acknowledgement:" << "<br />"
              << "'This product includes software developed by the OpenGeoSys Community.'</li>"
              << "<li>Neither the name of the OpenGeoSys Community nor the"
              << "names of its contributors may be used to endorse or promote products"
              << "derived from this software without specific prior written permission.</li>"
              << "<li>Attribute the OpenGeoSys Community, preferably citing an appropriate"
              << "paper listed on the OpenGeoSys Community homepage:"
              << "<a href='http://www.opengeosys.org/papers'>http://www.opengeosys.org/papers</a></li>"
           << "</ol></p>"
              << "<p>THIS SOFTWARE IS PROVIDED BY THE OpenGeoSys Community ''AS IS'' AND ANY"
              << "EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED"
              << "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE"
              << "DISCLAIMED. IN NO EVENT SHALL THE OpenGeoSys Community BE LIABLE FOR ANY"
              << "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES"
              << "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;"
              << "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND"
              << "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT"
              << "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS"
              << "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>";

    this->licenseTextBrowser->setHtml(*(stream.string()));
}

void LicenseDialog::on_okPushButton_pressed()
{
    this->done(QDialog::Accepted);
}

