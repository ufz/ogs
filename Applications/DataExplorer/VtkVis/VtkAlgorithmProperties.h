/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-23
 * \brief  Definition of the VtkAlgorithmProperties class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include <QList>
#include <QMap>
#include <QObject>
#include <QString>
#include <QVariant>

class vtkProperty;
class vtkTexture;
class vtkLookupTable;


#define ogsUserPropertyMacro(name,type) \
        virtual void Set ## name (type _arg) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " # name " to " << \
                _arg); \
        if (this->name != _arg)  \
        { \
            this->name = _arg; \
            this->Modified(); \
            (*(this->_algorithmUserProperties))[QString(# name)] = QVariant(_arg); \
        } \
    } \
\
        type name;
// End of ogsUserPropertyMacro

#define ogsUserVec2PropertyMacro(name,type) \
        virtual void Set ## name (type _arg1, type _arg2) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " << \
                # name " to (" << \
                _arg1 << "," << _arg2 << ")"); \
        if ((this->name[0] != _arg1) || (this->name[1] != _arg2)) \
        { \
            this->name[0] = _arg1; \
            this->name[1] = _arg2; \
            this->Modified(); \
            QList<QVariant> list; \
            list.push_back(QVariant(_arg1)); \
            list.push_back(QVariant(_arg2)); \
            (*(this->_algorithmUserVectorProperties))[QString(# name)] = list; \
        } \
    } \
\
        virtual void Set ## name (type _arg[2]) \
    { \
        this->Set ## name (_arg[0], _arg[1]); \
    } \
\
        type name[2];
// End of ogsUserVec2PropertyMacro

#define ogsUserVec3PropertyMacro(name,type) \
        virtual void Set ## name (type _arg1, type _arg2, type _arg3) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " << \
                # name " to (" << \
                _arg1 << "," << _arg2 << "," << _arg3 << ")"); \
        if ((this->name[0] != _arg1) || (this->name[1] != _arg2) || (this->name[2] != _arg3)) \
        { \
            this->name[0] = _arg1; \
            this->name[1] = _arg2; \
            this->name[2] = _arg3; \
            this->Modified(); \
            QList<QVariant> list; \
            list.push_back(QVariant(_arg1)); \
            list.push_back(QVariant(_arg2)); \
            list.push_back(QVariant(_arg3)); \
            (*(this->_algorithmUserVectorProperties))[QString(# name)] = list; \
        } \
    } \
\
        virtual void Set ## name (type _arg[3]) \
    { \
        this->Set ## name (_arg[0], _arg[1], _arg[2]); \
    } \
\
        type name[3];
// End of ogsUserVec3PropertyMacro

#define ogsUserVec4PropertyMacro(name,type) \
        virtual void Set ## name (type _arg1, type _arg2, type _arg3, type _arg4) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " << \
                # name " to (" << \
                _arg1 << "," << _arg2 << "," << _arg3 << "," << _arg4 << ")"); \
        if ((this->name[0] != _arg1) || (this->name[1] != _arg2) || \
            (this->name[2] != _arg3) || (this->name[3] != _arg4)) \
        { \
            this->name[0] = _arg1; \
            this->name[1] = _arg2; \
            this->name[2] = _arg3; \
            this->name[3] = _arg4; \
            this->Modified(); \
            QList<QVariant> list; \
            list.push_back(QVariant(_arg1)); \
            list.push_back(QVariant(_arg2)); \
            list.push_back(QVariant(_arg3)); \
            list.push_back(QVariant(_arg4)); \
            (*(this->_algorithmUserVectorProperties))[QString(# name)] = list; \
        } \
    } \
\
        virtual void Set ## name (type _arg[4]) \
    { \
        this->Set ## name (_arg[0], _arg[1], _arg[2], _arg[3]); \
    } \
\
        type name[4];
// End of ogsUserVec4PropertyMacro

/**
 * \brief Contains properties for the visualization of objects as VtkVisPipelineItems
 */
class VtkAlgorithmProperties : public QObject
{
    Q_OBJECT

public:
    /// Constructor (sets default values)
    VtkAlgorithmProperties(QObject* parent = nullptr);

    ~VtkAlgorithmProperties() override;

    /// @brief Returns the vtk properties
    vtkProperty* GetProperties() const { return _property; }

    /// @brief Returns a texture (if one has been assigned).
    vtkTexture* GetTexture() { return _texture; }

    /// @brief Sets a texture for the VtkVisPipelineItem.
    void SetTexture(vtkTexture* t) { _texture = t; }

    /// @brief Returns the colour lookup table (if one has been assigned).
    vtkLookupTable* GetLookupTable(const QString& array_name);

    /// @brief Removes the lookup table for the given scalar.
    void RemoveLookupTable(const QString& array_name);

    /// @brief Sets a colour lookup table for the given scalar array of the VtkVisPipelineItem.
    void SetLookUpTable(const QString &array_name, vtkLookupTable* lut);

    /// Loads a predefined color lookup table from a file for the specified scalar array.
    void SetLookUpTable(const QString &array_name, const QString &filename);

    /// @brief Returns the scalar visibility.
    bool GetScalarVisibility() const { return _scalarVisibility; }

    /// @brief Sets the scalar visibility.
    void SetScalarVisibility(bool on);

    /// @brief Returns the name. This is set to the file path if it is a source algorithm.
    QString GetName() const { return _name; }

    /// @brief Sets the name.
    void SetName(QString name) { _name = name; }

    /// @brief Is this algorithm removable from the pipeline (view).
    bool IsRemovable() const { return _removable; }

    /// @brief Returns a map of user properties.
    QMap<QString, QVariant>* GetAlgorithmUserProperties() const {
        return _algorithmUserProperties;
    }

    /// @brief Returns a map of vector user properties.
    QMap<QString, QList<QVariant> >* GetAlgorithmUserVectorProperties() const {
        return _algorithmUserVectorProperties;
    }

    /// @brief Sets a user property. This should be implemented by subclasses.
    virtual void SetUserProperty(QString name, QVariant value)
    {
        (*_algorithmUserProperties)[name] = value;
    }

    /// @brief Returns the value of a user property.
    QVariant GetUserProperty(QString name) const;

    /// @brief Sets a vector user property. This should be implemented by subclasses.
    virtual void SetUserVectorProperty(QString name, QList<QVariant> values)
    {
        (*_algorithmUserVectorProperties)[name] = values;
    }

    /// @brief Returns a list of values of a vector user property.
    QList<QVariant> GetUserVectorProperty(QString name) const;

    /// @brief Set the active attribute
    void SetActiveAttribute(QString name);

    /// @brief Returns the desired active attribute.
    QString GetActiveAttribute() const { return _activeAttributeName; }


protected:

    // Properties set on vtkActor
    vtkProperty* _property;
    vtkTexture* _texture;

    // Properties set on vtkMapper
    bool _scalarVisibility;
    std::map<QString, vtkLookupTable*> _lut;

    // Properties used in the GUI
    QString _name;
    QString _activeAttributeName;
    bool _removable;

    QMap<QString, QVariant>* _algorithmUserProperties;
    QMap<QString, QList<QVariant> >* _algorithmUserVectorProperties;

signals:
    void ScalarVisibilityChanged(bool on);
};
