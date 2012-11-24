/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file OptionGroup.h
 *
 * Created on 2012-09-25 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <map>
#include <ostream>
#include <vector>

#include "StringTools.h"
#include "OptionNode.h"
#include "OptionLeaf.h"

namespace BaseLib
{

/**
 * \brief Options represents a collection of key-value with hierarchical data structure.
 */
class OptionGroup : public OptionNode
{
private:
    typedef std::multimap<std::string, OptionNode*> DictionaryType;
    typedef std::pair<std::string, OptionNode*> PairType;
    typedef std::pair<DictionaryType::iterator, DictionaryType::iterator> pair_of_iterator;
    typedef std::pair<DictionaryType::const_iterator, DictionaryType::const_iterator> const_pair_of_iterator;
public:
    typedef DictionaryType::iterator iterator;
    typedef DictionaryType::const_iterator const_iterator;

    //_______________________________________________________________
    // Constructor
    //_______________________________________________________________
    /**
     * Constructor for a root node
     */
    OptionGroup() : _parent(NULL) {};

    /**
     * Create object with specifying a parent node
     * @param parent    a pointer to a parent node
     */
    explicit OptionGroup(OptionGroup* parent) : _parent(parent) {};

    /**
     *
     */
    virtual ~OptionGroup();


    //_______________________________________________________________
    /// add new option group 
    OptionGroup* addSubGroup(const std::string &key);

    /// return if a subgroup with the name exists
    bool hasSubGroup(const std::string &key) const;

    /// get option with the given key if it exists.
    OptionGroup* getSubGroup(const std::string &key);

    /// get option with the given key if it exists.
    const OptionGroup* getSubGroup(const std::string &key) const;

    /**
     * get parent option group
     * @return a pointer to parent
     */
    OptionGroup* getParentGroup() {return _parent;};


    //_______________________________________________________________
    /// add new option
    void addOption(const std::string &key, const std::string &v);

    /// add option as number
    template<typename T>
    inline void addOptionAsNum(const std::string &key, T v)
    {
        _dictionary.insert(PairType(key, new OptionLeaf<T>(v)));
    }

    /// check if there is a value with the given key
    bool hasOption(const std::string &key) const;

    /// get value as string
    const std::string getOption(const std::string &key) const;

    /// get raw data
    template<typename T>
    inline const T getOptionAsNum(const std::string &key) const
    {
        DictionaryType::const_iterator itr = _dictionary.find(key);
        if (itr==_dictionary.end() || !itr->second->isValue()) {
            return getDummy<T>();
        } else {
           if (itr->second->isString()) {
                return BaseLib::str2number<T>(static_cast<OptionLeaf<std::string>*>(itr->second)->getValue());
            } else {
                return static_cast<OptionLeaf<T>*>(itr->second)->getValue();
            }
        }
    }

    //_______________________________________________________________
    /// print out the content
    virtual void printout (std::ostream &os, size_t depth=0) const;

    //_______________________________________________________________
    /**
     * return iterator of the begin
     * @return
     */
    const_iterator begin() const { return _dictionary.begin(); };

    /**
     * return iterator of the end
     * @return
     */
    const_iterator end() const { return _dictionary.end(); };

    /// get a list of sub-group having the given key
    std::vector<const OptionGroup*> getSubGroupList(const std::string &key) const;

    /// get a list of values related to the given key
    template<typename T>
    std::vector<T> getOptionList(const std::string &key) const
    {
        // iterator for leaf
        const_pair_of_iterator _leaf_range = _dictionary.equal_range(key);
        DictionaryType::const_iterator itr;

        std::vector<T> list_op;

        for (itr = _leaf_range.first; itr!=_leaf_range.second; ++itr) {
            if (itr->second->isValue())
                list_op.push_back(static_cast<OptionLeaf<T>*>(itr->second)->getValue());
        }

        return list_op;
    }

private:
    // a dictionary storing pairs of key and OptionNodes
    DictionaryType _dictionary;
    // parent node
    OptionGroup* _parent;
};

}
