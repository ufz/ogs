/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-18
 * \brief  Implementation of the VtkVisTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkVisTabWidget.h"

#include "VtkAlgorithmPropertyCheckbox.h"
#include "VtkAlgorithmPropertyLineEdit.h"
#include "VtkAlgorithmPropertyVectorEdit.h"
#include "VtkColorByHeightFilter.h"
#include "VtkCompositeColorByHeightFilter.h"
#include "VtkVisImageItem.h"
#include "VtkVisPipelineItem.h"

#include <logog/include/logog.hpp>

#include <vtkActor.h>
#include <vtkImageChangeInformation.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

VtkVisTabWidget::VtkVisTabWidget( QWidget* parent /*= 0*/ )
    : QWidget(parent), _item(nullptr)
{
    setupUi(this);

    this->scaleZ->setValidator(new QDoubleValidator(0, 100, 8, this));

    this->transX->setValidator(new QDoubleValidator(this));
    this->transY->setValidator(new QDoubleValidator(this));
    this->transZ->setValidator(new QDoubleValidator(this));

    connect(this->vtkVisPipelineView, SIGNAL(requestViewUpdate()),
            this, SIGNAL(requestViewUpdate()));

    connect(this->vtkVisPipelineView, SIGNAL(itemSelected(VtkVisPipelineItem*)),
            this, SLOT(setActiveItem(VtkVisPipelineItem*)));

    connect(this->activeScalarComboBox, SIGNAL(currentIndexChanged(const QString &)),
            this, SLOT(SetActiveAttributeOnItem(const QString &)));
}

void VtkVisTabWidget::on_arrayResetPushButton_clicked()
{
    VtkAlgorithmProperties* props = _item->getVtkProperties();
    const QString selected_array_name = this->activeScalarComboBox->currentText();
    props->RemoveLookupTable(selected_array_name);
    _item->SetActiveAttribute(selected_array_name);
}

void VtkVisTabWidget::setActiveItem( VtkVisPipelineItem* item )
{
    if (item)
    {
        _item = item;
        transformTabWidget->setEnabled(true);

        vtkTransformFilter* transform_filter = dynamic_cast<vtkTransformFilter*>(_item->transformFilter());
        if (transform_filter) // if data set
        {
            actorPropertiesGroupBox->setEnabled(true);
            vtkProperty* vtkProps = static_cast<vtkActor*>(_item->actor())->GetProperty();
            diffuseColorPickerButton->setColor(vtkProps->GetDiffuseColor());
            visibleEdgesCheckBox->setChecked(vtkProps->GetEdgeVisibility());
            edgeColorPickerButton->setColor(vtkProps->GetEdgeColor());
            opacitySlider->setValue((int)(vtkProps->GetOpacity() * 100.0));

            vtkTransform* transform =
                    static_cast<vtkTransform*>(transform_filter->GetTransform());
            if (transform)
            {
                double scale[3];
                transform->GetScale(scale);
                double trans[3];
                transform->GetPosition(trans);

                //switch signals off for just filling in text-boxes after clicking on an item
                this->scaleZ->blockSignals(true);
                this->transX->blockSignals(true);
                this->transY->blockSignals(true);
                this->transZ->blockSignals(true);
                this->scaleZ->setText(QString::number(scale[2]));
                this->transX->setText(QString::number(trans[0] / scale[0]));
                this->transY->setText(QString::number(trans[1] / scale[1]));
                this->transZ->setText(QString::number(trans[2] / scale[2]));
                this->scaleZ->blockSignals(false);
                this->transX->blockSignals(false);
                this->transY->blockSignals(false);
                this->transZ->blockSignals(false);
                //switch signals back on
            }
            this->buildScalarArrayComboBox(_item);

            // Set to last active attribute
            QString activeAttribute = _item->GetActiveAttribute();
            if (activeAttribute.length() > 0)
                for (int i = 0; i < this->activeScalarComboBox->count(); i++)
                {
                    QString itemText = this->activeScalarComboBox->itemText(i);
                    if (itemText.compare(activeAttribute) == 0)
                    {
                        this->activeScalarComboBox->setCurrentIndex(i);
                        break;
                    }
                }
        }
        else // if image
        {
            const VtkVisImageItem* img = static_cast<VtkVisImageItem*>(_item);
            actorPropertiesGroupBox->setEnabled(false);
            vtkImageChangeInformation* transform = static_cast<vtkImageChangeInformation*>(img->transformFilter());
            double trans[3];
            transform->GetOriginTranslation(trans);
            this->transX->blockSignals(true);
            this->transY->blockSignals(true);
            this->transZ->blockSignals(true);
            this->transX->setText(QString::number(trans[0]));
            this->transY->setText(QString::number(trans[1]));
            this->transZ->setText(QString::number(trans[2]));
            this->transX->blockSignals(false);
            this->transY->blockSignals(false);
            this->transZ->blockSignals(false);
        }

        this->buildProportiesDialog(item);

        emit requestViewUpdate();
    }
    else
    {
        actorPropertiesGroupBox->setEnabled(false);
        transformTabWidget->setEnabled(false);
        this->activeScalarComboBox->clear();
    }
}

void VtkVisTabWidget::on_diffuseColorPickerButton_colorPicked( QColor color )
{
    static_cast<vtkActor*>(_item->actor())->GetProperty()->SetDiffuseColor(
            color.redF(), color.greenF(), color.blueF());

    emit requestViewUpdate();
}

void VtkVisTabWidget::on_visibleEdgesCheckBox_stateChanged( int state )
{
    if (state == Qt::Checked)
    {
        static_cast<vtkActor*>(_item->actor())->GetProperty()->SetEdgeVisibility(1);
        edgeColorPickerButton->setEnabled(true);
    }
    else
    {
        static_cast<vtkActor*>(_item->actor())->GetProperty()->SetEdgeVisibility(0);
        edgeColorPickerButton->setEnabled(false);
    }

    emit requestViewUpdate();
}

void VtkVisTabWidget::on_edgeColorPickerButton_colorPicked( QColor color )
{
    static_cast<vtkActor*>(_item->actor())->GetProperty()->SetEdgeColor(
            color.redF(), color.greenF(), color.blueF());
    emit requestViewUpdate();
}

void VtkVisTabWidget::on_opacitySlider_sliderMoved( int value )
{
    static_cast<vtkActor*>(_item->actor())->GetProperty()->SetOpacity(value / 100.0);
    emit requestViewUpdate();
}

void VtkVisTabWidget::on_scaleZ_textChanged(const QString &text)
{
    bool ok = true;
    double scale = text.toDouble(&ok);

    // If z scale becomes zero, the object becomes invisible
    if (ok && scale != 0.0)
    {
        _item->setScale(1.0, 1.0, scale);

        for (int i = 0; i < _item->childCount(); i++)
        {
            VtkVisPipelineItem* childItem = _item->child(i);
            if (childItem)
            {
                VtkCompositeColorByHeightFilter* colorFilter =
                        dynamic_cast<VtkCompositeColorByHeightFilter*>
                        (childItem->compositeFilter());
                if (colorFilter)
                    VtkColorByHeightFilter::SafeDownCast(
                            colorFilter->GetOutputAlgorithm())->
                    SetTableRangeScaling(scale);
            }
        }

        emit requestViewUpdate();
    }
}

void VtkVisTabWidget::translateItem()
{
    bool okX(true), okY(true), okZ(true);
    double trans[3];

    trans[0] = transX->text().toDouble(&okX);
    trans[1] = transY->text().toDouble(&okY);
    trans[2] = transZ->text().toDouble(&okZ);

    if (okX && okY && okZ)
    {
        _item->setTranslation(trans[0], trans[1], trans[2]);
        emit requestViewUpdate();
    }
}

void VtkVisTabWidget::buildProportiesDialog(VtkVisPipelineItem* item)
{
    QFormLayout* layout = static_cast<QFormLayout*>(this->scrollAreaWidgetContents->layout());
    while(layout->count())
        delete layout->takeAt(0)->widget();

    QMap<QString, QVariant>* propMap = nullptr;
    QMap<QString, QList<QVariant>>* propVecMap = nullptr;
    VtkAlgorithmProperties* algProps = item->getVtkProperties();

    if (algProps == nullptr)
        WARN("VtkAlgorithmProperties null!")

    // Retrieve algorithm properties
    if (item->compositeFilter())
    {
        propMap = item->compositeFilter()->GetAlgorithmUserProperties();
        propVecMap = item->compositeFilter()->GetAlgorithmUserVectorProperties();
    }
    else
    {
        if (algProps)
        {
            propMap = algProps->GetAlgorithmUserProperties();
            propVecMap = algProps->GetAlgorithmUserVectorProperties();
        }
    }

    // Select appropriate GUI element and set connect for each property
    if (propMap)
    {
        QMapIterator<QString, QVariant> i(*propMap);
        while (i.hasNext())
        {
            i.next();
            QString key = i.key();
            QVariant value = i.value();

            VtkAlgorithmPropertyLineEdit* lineEdit;
            VtkAlgorithmPropertyCheckbox* checkbox;
            switch (value.type())
            {
            case QVariant::Double:
                lineEdit =
                        new VtkAlgorithmPropertyLineEdit(QString::number(
                                                                 value.toDouble()),
                                                         key, QVariant::Double,
                                                         algProps);
                connect(lineEdit, SIGNAL(editingFinished()), this,
                        SIGNAL(requestViewUpdate()));
                layout->addRow(key, lineEdit);
                break;

            case QVariant::Int:
                lineEdit =
                        new VtkAlgorithmPropertyLineEdit(QString::number(
                                                                 value.toInt()),
                                                         key, QVariant::Int,
                                                         algProps);
                connect(lineEdit, SIGNAL(editingFinished()), this,
                        SIGNAL(requestViewUpdate()));
                layout->addRow(key, lineEdit);
                break;

            case QVariant::Bool:
                checkbox = new VtkAlgorithmPropertyCheckbox(
                        value.toBool(), key, algProps);
                connect(checkbox, SIGNAL(stateChanged(int)), this,
                        SIGNAL(requestViewUpdate()));
                layout->addRow(key, checkbox);
                break;

            default:
                break;
            }
        }
    }

    if (propVecMap)
    {
        QMapIterator<QString, QList<QVariant> > i(*propVecMap);
        while (i.hasNext())
        {
            i.next();
            QString key = i.key();
            QList<QVariant> values = i.value();

            VtkAlgorithmPropertyVectorEdit* vectorEdit;
            if (values.size() > 0)
            {
                QList<QString> valuesAsString;
                foreach (QVariant variant, values)
                valuesAsString.push_back(variant.toString());

                vectorEdit = new VtkAlgorithmPropertyVectorEdit(valuesAsString, key,
                                                                values.front().type(),
                                                                algProps);
                connect(vectorEdit, SIGNAL(editingFinished()), this,
                        SIGNAL(requestViewUpdate()));
                layout->addRow(key, vectorEdit);
            }
        }
    }
}

void VtkVisTabWidget::buildScalarArrayComboBox(VtkVisPipelineItem* item)
{
    QStringList dataSetAttributesList = item->getScalarArrayNames();
    dataSetAttributesList.push_back("Solid Color"); // all scalars switched off
    this->activeScalarComboBox->blockSignals(true);
    this->activeScalarComboBox->clear();
    this->activeScalarComboBox->insertItems(0, dataSetAttributesList);
    this->activeScalarComboBox->blockSignals(false);
    QString active_array_name = item->GetActiveAttribute();
    QList<QString>::iterator it = dataSetAttributesList.begin();
    if (active_array_name.length() == 0)
        item->SetActiveAttribute(*it);
    else
    {
        int idx(0);
        for (it=dataSetAttributesList.begin(); it!=dataSetAttributesList.end(); ++it)
            if (active_array_name.compare((*it).right((*it).length()-2))==0)
            {
                this->activeScalarComboBox->setCurrentIndex(idx);
                break;
            }
            else idx++;
    }
}

void VtkVisTabWidget::SetActiveAttributeOnItem( const QString &name )
{
    _item->SetActiveAttribute(name);
    emit requestViewUpdate();
}

