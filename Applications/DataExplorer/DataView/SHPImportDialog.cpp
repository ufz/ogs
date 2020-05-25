/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-25
 * \brief  Implementation of the SHPImportDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoLib/GEOObjects.h"
#include "OGSError.h"
#include "SHPImportDialog.h"
#include "SHPInterface.h"

#include <QDialogButtonBox>
#include <QFileInfo>
#include <QLabel>
#include <QLineEdit>
#include <QRadioButton>
#include <QVBoxLayout>

SHPImportDialog::SHPImportDialog(std::string filename,
                                 GeoLib::GEOObjects& geo_objects,
                                 std::string const& gmsh_path,
                                 QDialog* parent)
    : QDialog(parent),
      buttonBox_(nullptr),
      layout_(nullptr),
      shpContentLabel_(nullptr),
      nameLabel_(nullptr),
      listName_(new QLineEdit()),
      choice1_(nullptr),
      choice2_(nullptr),
      filename_(std::move(filename)),
      fileType_(0),
      shpInterface_(new FileIO::SHPInterface(geo_objects)),
      gmsh_path_(gmsh_path)
{
    setupDialog();
    show();
}

SHPImportDialog::~SHPImportDialog()
{
    delete shpInterface_;
    delete buttonBox_;
    delete listName_;
    delete choice1_;
    delete choice2_;
    delete shpContentLabel_;
    delete nameLabel_;
    delete layout_;
}

void SHPImportDialog::setupDialog()
{
    layout_ = new QGridLayout(this);
    int shpType = 0;
    int numberOfEntities = 0;
    QString type = "";

    setWindowTitle("Import SHP File");

    if (shpInterface_->readSHPInfo(filename_, shpType, numberOfEntities))
    {
        if ((shpType - 1) % 10 == 0)
        {
            type = "points";
        }
        if ((shpType - 3) % 10 == 0)
        {
            type = "polylines";
        }
        if ((shpType - 5) % 10 == 0)
        {
            type = "polygons";
        }
        if ((shpType - 8) % 10 == 0)
        {
            type = "multipoints";
        }
        if (shpType == 31)
        {
            type = "TIN elements";
        }

        shpContentLabel_ =
                new QLabel("The selected file contains " + QString::number(
                                   numberOfEntities) + " " + type, this);
        nameLabel_ = new QLabel("List Name: ", this);

        QFileInfo fi(QString::fromStdString(filename_));
        listName_->setText(fi.baseName());

        if ((shpType - 1) % 10 == 0 && shpType != 31) // Points
        {
            choice1_ = new QRadioButton("Read as Geometry Points");
            choice2_ = new QRadioButton("Read as Station Points");
            choice1_->toggle(); // default choice
            layout_->addWidget( shpContentLabel_ );
            layout_->addWidget( choice1_ );
            layout_->addWidget( choice2_ );
            layout_->addWidget( nameLabel_ );
            layout_->addWidget( listName_ );
            fileType_ = 1;
        }
        else if ((shpType - 3) % 10 == 0 || (shpType - 5) % 10 == 0) // Polylines
        {
            choice1_ = new QRadioButton("Read Polylines only");
            choice2_ = new QRadioButton("Read Polylines/Surfaces");
            if ((shpType - 3) % 10 == 0)
            {
                choice2_->setDisabled(true);  // disable polygon-choice if file
                                              // contains only polylines
            }
            choice1_->toggle(); // default choice
            layout_->addWidget( shpContentLabel_ );
            layout_->addWidget( choice1_ );
            layout_->addWidget( choice2_ );
            layout_->addWidget( nameLabel_ );
            layout_->addWidget( listName_ );
            fileType_ = 2;
        }
        else
        {
            nameLabel_->setText("This element type is currently not supported.");
            layout_->addWidget( shpContentLabel_ );
            layout_->addWidget( nameLabel_ );
        }

        buttonBox_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
        connect(buttonBox_, SIGNAL(accepted()), this, SLOT(accept()));
        connect(buttonBox_, SIGNAL(rejected()), this, SLOT(reject()));
        layout_->addWidget( buttonBox_ );

        setLayout(layout_);
    }
    else
    {
        OGSError::box("Error reading shapefile!");
    }
}

void SHPImportDialog::accept()
{
    QString list_name(listName_->text());
    if (list_name.compare("") == 0)
    {
        OGSError::box("Please insert a name for the data in this file.");
        return;
    }

    if (fileType_ == 1 && choice1_->isChecked())
    {
        shpInterface_->readSHPFile(filename_,
                                   FileIO::SHPInterface::OGSType::POINT,
                                   list_name.toStdString(), gmsh_path_);
    }
    if (fileType_ == 1 && choice2_->isChecked())
    {
        shpInterface_->readSHPFile(filename_,
                                   FileIO::SHPInterface::OGSType::STATION,
                                   list_name.toStdString(), gmsh_path_);
    }
    if (fileType_ == 2 && choice1_->isChecked())
    {
        shpInterface_->readSHPFile(filename_,
                                   FileIO::SHPInterface::OGSType::POLYLINE,
                                   list_name.toStdString(), gmsh_path_);
    }
    if (fileType_ == 2 && choice2_->isChecked())
    {
        shpInterface_->readSHPFile(filename_,
                                   FileIO::SHPInterface::OGSType::POLYGON,
                                   list_name.toStdString(), gmsh_path_);
    }
    emit shpLoaded(list_name);

    this->done(QDialog::Accepted);
}

void SHPImportDialog::reject()
{
    this->done(QDialog::Rejected);
}
