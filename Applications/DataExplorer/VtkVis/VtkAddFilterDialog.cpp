/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-23
 * \brief  Implementation of the VtkAddFilterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkAddFilterDialog.h"

#include <logog/include/logog.hpp>

#include "VtkCompositeFilter.h"
#include "VtkFilterFactory.h"
#include "VtkVisImageItem.h"
#include "VtkVisPipeline.h"
#include "VtkVisPipelineItem.h"
#include "VtkVisPointSetItem.h"

#include <vtkContourFilter.h>
#include <vtkOutlineFilter.h>
#include <vtkTransformFilter.h>

#include <QModelIndex>

VtkAddFilterDialog::VtkAddFilterDialog( VtkVisPipeline &pipeline,
                                        QModelIndex parentIndex,
                                        QDialog* parent /*= 0*/ )
    : QDialog(parent), _pipeline(pipeline), _parentIndex(parentIndex)
{
    setupUi(this);
    filterListWidget->setSelectionMode(QAbstractItemView::SingleSelection);

    auto* parentItem =
        static_cast<VtkVisPipelineItem*>(_pipeline.getItem(parentIndex));
    vtkDataObject* parentDataObject = parentItem->algorithm()->GetOutputDataObject(0);
    int parentDataObjectType = parentDataObject->GetDataObjectType();

    QVector<VtkFilterInfo> filterList = VtkFilterFactory::GetFilterList();
    foreach(VtkFilterInfo filter, filterList)
    {
        // Check for suitable filters (see vtkDataSet inheritance diagram)
        int inputType = filter.inputDataObjectType;
        if ((inputType == parentDataObjectType) ||
            (inputType == VTK_POINT_SET && parentDataObjectType != VTK_IMAGE_DATA) ||
            (inputType == VTK_IMAGE_DATA &&
             (parentDataObjectType == VTK_STRUCTURED_POINTS || parentDataObjectType ==
              VTK_UNIFORM_GRID)))

            new QListWidgetItem(filter.readableName, filterListWidget);
    }

    // On double clicking an item the dialog gets accepted
    connect(filterListWidget,SIGNAL(itemDoubleClicked(QListWidgetItem*)),
            this->buttonBox,SIGNAL(accepted()));
}

void VtkAddFilterDialog::on_buttonBox_accepted()
{
    QVector<VtkFilterInfo> filterList = VtkFilterFactory::GetFilterList();
    QString filterName;
    foreach(VtkFilterInfo filter, filterList)
    {
        if (filter.readableName.compare(filterListWidget->currentItem()->text()) == 0)
        {
            filterName = filter.name;
            break;
        }
    }
    auto* parentItem =
        static_cast<VtkVisPipelineItem*>(_pipeline.getItem(_parentIndex));
    QList<QVariant> itemData;
    itemData << filterListWidget->currentItem()->text() << true;

    VtkCompositeFilter* filter;
    if (dynamic_cast<VtkVisImageItem*>(parentItem))
        filter = VtkFilterFactory::CreateCompositeFilter(filterName, parentItem->algorithm());
    else
        filter = VtkFilterFactory::CreateCompositeFilter(filterName,
                                                         parentItem->transformFilter());

    VtkVisPipelineItem* item;
    if (filter)
    {
        if (filter->GetOutputDataObjectType() == VTK_IMAGE_DATA)
            item = new VtkVisImageItem(filter, parentItem, itemData);
        else
            item = new VtkVisPointSetItem(filter, parentItem, itemData);
    }
    else
    {
        vtkAlgorithm* algorithm = VtkFilterFactory::CreateSimpleFilter(filterName);
        if (algorithm)
            item = new VtkVisPointSetItem(algorithm, parentItem, itemData);
        else
        {
            ERR("VtkFilterFactory cannot create %s.", filterName.toStdString().c_str());
            return;
        }
    }
    _pipeline.addPipelineItem(item, _parentIndex);
}

void VtkAddFilterDialog::on_filterListWidget_currentRowChanged( int currentRow )
{
    foreach(VtkFilterInfo filter, VtkFilterFactory::GetFilterList())
    {
        if (filter.readableName.compare(filterListWidget->item(currentRow)->text()) == 0)
        {
            QString desc = filter.description;
            desc = desc + QString("\n\nThis filter outputs ") +
                   filter.OutputDataObjectTypeAsString() +
                   QString("\n\nFilter class name: ") +
                   filter.name;

            this->filterDescTextEdit->setText(desc);
            continue;
        }
    }
}
