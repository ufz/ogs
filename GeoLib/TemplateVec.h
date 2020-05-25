/**
 * @file TemplateVec.h
 * @date 2010-02-26
 * @author Thomas Fischer
 * @brief Definition of the GeoLib::TemplateVec class.
 *
 * @copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "BaseLib/Logging.h"

#include "BaseLib/Error.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 *
 * \brief The class TemplateVec takes a unique name and manages
 * a std::vector of pointers to data elements of type T.
 *
 * Instances are classes PolylineVec and SurfaceVec.
 * */
template <class T> class TemplateVec
{
protected:
    using NameIdPair = std::pair<std::string, std::size_t>;
    using NameIdMap = std::map<std::string, std::size_t>;

public:
    /**
     * Constructor of class TemlateVec.
     * @param name unique name of the project the elements belonging to.
     * In order to access the data elements a unique name is required.
     * @param data_vec Vector of data elements.
     * @attention{TemplateVec will take the ownership of the vector
     * and also its elements,
     * i.e. delete its elements and delete the vector itself!}
     * @param elem_name_map Names of data elements can be given by a
     * std::map<std::string, std::size_t>. Here the std::string is the name
     * of the element and the value for std::size_t stands for an index in
     * the data_vec.
     */
    TemplateVec(std::string name, std::unique_ptr<std::vector<T*>> data_vec,
                std::unique_ptr<NameIdMap> elem_name_map = nullptr)
        : name_(std::move(name)),
          data_vec_(std::move(data_vec)),
          name_id_map_(std::move(elem_name_map))
    {
        if (data_vec_ == nullptr)
        {
            OGS_FATAL("Constructor TemplateVec: vector of data elements is a nullptr.");
        }

        if (!name_id_map_)
        {
            name_id_map_ = std::make_unique<NameIdMap>();
        }
    }

    /**
     * destructor, deletes all data elements
     */
    virtual ~TemplateVec ()
    {
        for (std::size_t k(0); k < size(); k++)
        {
            delete (*data_vec_)[k];
        }
    }

    /** sets the name of the vector of geometric objects
     * the data elements belonging to
     * \param n the name as standard string */
    void setName (const std::string & n) { name_ = n; }
    /**
     * the name, the data element belonging to
     * @return the name of the object
     */
    std::string getName () const { return name_; }

    /// Returns the begin of the name id mapping structure
    NameIdMap::const_iterator getNameIDMapBegin() const { return name_id_map_->cbegin(); }

    /// Returns the end of the name id mapping structure
    NameIdMap::const_iterator getNameIDMapEnd() const { return name_id_map_->cend(); }

    /**
     * @return the number of data elements
     */
    std::size_t size () const { return data_vec_->size(); }

    /**
     * get a pointer to a standard vector containing the data elements
     * @return the data elements
     */
    const std::vector<T*>* getVector () const { return data_vec_.get(); }

    /**
     * search the vector of names for the ID of the geometric element with the given name
     * @param name the name of the geometric element
     * @param id the id of the geometric element
     */
    bool getElementIDByName (const std::string& name, std::size_t &id) const
    {
        auto it (name_id_map_->find (name));

        if (it != name_id_map_->end())
        {
            id = it->second;
            return true;
        }
        return false;
    }

    /// Returns an element with the given name.
    const T* getElementByName (const std::string& name) const
    {
        std::size_t id;
        bool ret (getElementIDByName (name, id));
        if (ret)
        {
            return (*data_vec_)[id];
        }

        return nullptr;
    }

    /**
     * The method returns true if there is a name associated
     * with the given id, else method returns false.
     * @param id the id
     * @param element_name if a name associated with the id
     * is found name is assigned to element_name
     * @return if there is name associated with the given id:
     * true, else false
     */
    bool getNameOfElementByID (std::size_t id, std::string& element_name) const
    {
        // search in map for id
        auto it = findFirstElementByID(id);
        if (it == name_id_map_->end()) {
            return false;
        }
        element_name = it->first;
        return true;
    }

    /// Return the name of an element based on its ID.
    void setNameOfElementByID (std::size_t id, std::string const& element_name)
    {
        name_id_map_->insert(NameIdPair(element_name, id));
    }

    /**
     * The method returns true if the given element of type T
     * can be found and the element has a name, else method returns false.
     * @param data the data element, one wants to know the name
     * @param name the name of the data element (if the data element is
     * found and the data element has a name)
     * @return if element is found and has a name: true, else false
     */
    bool getNameOfElement (const T* data, std::string& name) const
    {
        for (std::size_t k(0); k < data_vec_->size(); k++)
        {
            if ((*data_vec_)[k] == data)
            {
                return getNameOfElementByID(k, name);
            }
        }

        return false;
    }

    /// Adds a new element to the vector.
    virtual void push_back (T* data_element, std::string const* const name = nullptr)
    {
        data_vec_->push_back (data_element);
        if (!name || name->empty())
        {
            return;
        }

        std::map<std::string, std::size_t>::const_iterator it(
            name_id_map_->find(*name)
        );
        if (it == name_id_map_->end()) {
            name_id_map_->insert(NameIdPair(*name, data_vec_->size() - 1));
        } else {
            WARN(
                "Name '{:s}' exists already. The object will be inserted "
                "without a name",
                name->c_str());
        }
    }

    /// Sets the given name for the element of the given ID.
    virtual void setNameForElement(std::size_t id, std::string const& name)
    {
        // Erase id if found in map.
        auto it = findFirstElementByID(id);
        if (it != name_id_map_->end())
        {
            name_id_map_->erase(it);
        }

        if (!name.empty()) {
            //insert new or revised name
            name_id_map_->insert(NameIdPair(name, id));
        }
    }

private:

    NameIdMap::const_iterator
    findFirstElementByID(std::size_t const& id) const
    {
        return std::find_if(name_id_map_->begin(), name_id_map_->end(),
            [id](NameIdPair const& elem) { return elem.second == id; });
    }

protected:
    /** copy constructor doesn't have an implementation */
    // compiler does not create a (possible unwanted) copy constructor
    TemplateVec (const TemplateVec &);
    /** assignment operator doesn't have an implementation */
    // this way the compiler does not create a (possible unwanted) assignment operator
    TemplateVec& operator= (const TemplateVec& rhs);

    /** the name of the object */
    std::string name_;

    /**
     * pointer to a vector of data elements
     */
    std::unique_ptr<std::vector <T*>> data_vec_;
    /**
     * store names associated with the element ids
     */
    std::unique_ptr<NameIdMap> name_id_map_;
};
} // end namespace GeoLib
