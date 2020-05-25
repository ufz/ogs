/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-23
 * \brief  Definition of the VtkAlgorithmProperties class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
        virtual void Set ## name (type arg_) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " # name " to " << \
                arg_); \
        if (this->name != arg_)  \
        { \
            this->name = arg_; \
            this->Modified(); \
            (*(this->algorithmUserProperties_))[QString(# name)] = QVariant(arg_); \
        } \
    } \
\
        type name;
// End of ogsUserPropertyMacro

#define ogsUserVec2PropertyMacro(name,type) \
        virtual void Set ## name (type arg1_, type arg2_) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " << \
                # name " to (" << \
                arg1_ << "," << arg2_ << ")"); \
        if ((this->name[0] != arg1_) || (this->name[1] != arg2_)) \
        { \
            this->name[0] = arg1_; \
            this->name[1] = arg2_; \
            this->Modified(); \
            QList<QVariant> list; \
            list.push_back(QVariant(arg1_)); \
            list.push_back(QVariant(arg2_)); \
            (*(this->algorithmUserVectorProperties_))[QString(# name)] = list; \
        } \
    } \
\
        virtual void Set ## name (type arg_[2]) \
    { \
        this->Set ## name (arg_[0], arg_[1]); \
    } \
\
        type name[2];
// End of ogsUserVec2PropertyMacro

#define ogsUserVec3PropertyMacro(name,type) \
        virtual void Set ## name (type arg1_, type arg2_, type arg3_) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " << \
                # name " to (" << \
                arg1_ << "," << arg2_ << "," << arg3_ << ")"); \
        if ((this->name[0] != arg1_) || (this->name[1] != arg2_) || (this->name[2] != arg3_)) \
        { \
            this->name[0] = arg1_; \
            this->name[1] = arg2_; \
            this->name[2] = arg3_; \
            this->Modified(); \
            QList<QVariant> list; \
            list.push_back(QVariant(arg1_)); \
            list.push_back(QVariant(arg2_)); \
            list.push_back(QVariant(arg3_)); \
            (*(this->algorithmUserVectorProperties_))[QString(# name)] = list; \
        } \
    } \
\
        virtual void Set ## name (type arg_[3]) \
    { \
        this->Set ## name (arg_[0], arg_[1], arg_[2]); \
    } \
\
        type name[3];
// End of ogsUserVec3PropertyMacro

#define ogsUserVec4PropertyMacro(name,type) \
        virtual void Set ## name (type arg1_, type arg2_, type arg3_, type arg4_) \
    { \
        vtkDebugMacro( \
                << this->GetClassName() << " (" << this << "): setting " << \
                # name " to (" << \
                arg1_ << "," << arg2_ << "," << arg3_ << "," << arg4_ << ")"); \
        if ((this->name[0] != arg1_) || (this->name[1] != arg2_) || \
            (this->name[2] != arg3_) || (this->name[3] != arg4_)) \
        { \
            this->name[0] = arg1_; \
            this->name[1] = arg2_; \
            this->name[2] = arg3_; \
            this->name[3] = arg4_; \
            this->Modified(); \
            QList<QVariant> list; \
            list.push_back(QVariant(arg1_)); \
            list.push_back(QVariant(arg2_)); \
            list.push_back(QVariant(arg3_)); \
            list.push_back(QVariant(arg4_)); \
            (*(this->algorithmUserVectorProperties_))[QString(# name)] = list; \
        } \
    } \
\
        virtual void Set ## name (type arg_[4]) \
    { \
        this->Set ## name (arg_[0], arg_[1], arg_[2], arg_[3]); \
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
    explicit VtkAlgorithmProperties(QObject* parent = nullptr);

    ~VtkAlgorithmProperties() override;

    /// @brief Returns the vtk properties
    vtkProperty* GetProperties() const { return property_; }

    /// @brief Returns a texture (if one has been assigned).
    vtkTexture* GetTexture() { return texture_; }

    /// @brief Sets a texture for the VtkVisPipelineItem.
    void SetTexture(vtkTexture* t) { texture_ = t; }

    /// @brief Returns the colour lookup table (if one has been assigned).
    vtkLookupTable* GetLookupTable(const QString& array_name);

    /// @brief Removes the lookup table for the given scalar.
    void RemoveLookupTable(const QString& array_name);

    /// @brief Sets a colour lookup table for the given scalar array of the VtkVisPipelineItem.
    void SetLookUpTable(const QString &array_name, vtkLookupTable* lut);

    /// Loads a predefined color lookup table from a file for the specified scalar array.
    void SetLookUpTable(const QString &array_name, const QString &filename);

    /// @brief Returns the scalar visibility.
    bool GetScalarVisibility() const { return scalarVisibility_; }

    /// @brief Sets the scalar visibility.
    void SetScalarVisibility(bool on);

    /// @brief Returns the name. This is set to the file path if it is a source algorithm.
    QString GetName() const { return name_; }

    /// @brief Sets the name.
    void SetName(QString name) { name_ = name; }

    /// @brief Is this algorithm removable from the pipeline (view).
    bool IsRemovable() const { return removable_; }

    /// @brief Returns a map of user properties.
    QMap<QString, QVariant>* GetAlgorithmUserProperties() const {
        return algorithmUserProperties_;
    }

    /// @brief Returns a map of vector user properties.
    QMap<QString, QList<QVariant> >* GetAlgorithmUserVectorProperties() const {
        return algorithmUserVectorProperties_;
    }

    /// @brief Sets a user property. This should be implemented by subclasses.
    virtual void SetUserProperty(QString name, QVariant value)
    {
        (*algorithmUserProperties_)[name] = value;
    }

    /// @brief Returns the value of a user property.
    QVariant GetUserProperty(QString name) const;

    /// @brief Sets a vector user property. This should be implemented by subclasses.
    virtual void SetUserVectorProperty(QString name, QList<QVariant> values)
    {
        (*algorithmUserVectorProperties_)[name] = values;
    }

    /// @brief Returns a list of values of a vector user property.
    QList<QVariant> GetUserVectorProperty(QString name) const;

    /// @brief Set the active attribute
    void SetActiveAttribute(QString name);

    /// @brief Returns the desired active attribute.
    QString GetActiveAttribute() const { return activeAttributeName_; }


protected:

    // Properties set on vtkActor
    vtkProperty* property_;
    vtkTexture* texture_;

    // Properties set on vtkMapper
    bool scalarVisibility_;
    std::map<QString, vtkLookupTable*> lut_;

    // Properties used in the GUI
    QString name_;
    QString activeAttributeName_;
    bool removable_;

    QMap<QString, QVariant>* algorithmUserProperties_;
    QMap<QString, QList<QVariant> >* algorithmUserVectorProperties_;

signals:
    void ScalarVisibilityChanged(bool on);
};
