/**
 * \file
 * \author Lars Bilke
 * \date   2010-11-12
 * \brief  Definition of the QVtkDataSetMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QVTKDATASETMAPPER_H
#define QVTKDATASETMAPPER_H

#include <QObject>
#include <vtkDataSetMapper.h>

/// @brief Simply wraps vtkDataSetMapper as a Qt object to enable slot connections.
class QVtkDataSetMapper : public QObject, public vtkDataSetMapper
{
    Q_OBJECT

public:
    /// @brief Create new objects with New() because of VTKs reference counting.
    static QVtkDataSetMapper* New();

    vtkTypeMacro(QVtkDataSetMapper, vtkDataSetMapper);

    /// @brief Prints information about itself.
    void PrintSelf(ostream& os, vtkIndent indent) override;

public slots:
    /// @brief Sets the scalar visibility on this mapper.
    virtual void SetScalarVisibility(bool on);
    virtual void SetScalarVisibility(int on) { SetScalarVisibility(static_cast<bool>(on)); }

protected:
    /// @brief Constructor.
    QVtkDataSetMapper();

    /// @brief Destructor.
    virtual ~QVtkDataSetMapper();

private:
    QVtkDataSetMapper(const QVtkDataSetMapper&); // Not implemented.
    void operator=(const QVtkDataSetMapper&); // Not implemented
};

#endif // QVTKDATASETMAPPER_H
