/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

template <class T>
class InstanceCounter
{
public:
    InstanceCounter() {
        ++_num_constructed;
    }
    InstanceCounter(InstanceCounter<T> const&) {
        ++_num_copied;
    }
    InstanceCounter(InstanceCounter<T>&&) {
        ++_num_moved;
    }
    virtual ~InstanceCounter() {
        ++_num_destroyed;
    }

    static int getNumberOfConstructions() { return _num_constructed; }
    static int getNumberOfCopies() { return _num_copied; }
    static int getNumberOfMoves() { return _num_moved; }
    static int getNumberOfDestructions() { return _num_destroyed; }
    static int getNumberOfInstances()
    {
        return _num_constructed + _num_moved + _num_copied - _num_destroyed;
    }

    static void update(int& num_const, int& num_move, int& num_copy, int& num_inst)
    {
        num_const = getNumberOfConstructions();
        num_move = getNumberOfMoves();
        num_copy = getNumberOfCopies();
        num_inst = getNumberOfInstances();
    }

private:
    static int _num_constructed;
    static int _num_copied;
    static int _num_moved;
    static int _num_destroyed;
};

template <class T>
int InstanceCounter<T>::_num_constructed = 0;
template <class T>
int InstanceCounter<T>::_num_copied = 0;
template <class T>
int InstanceCounter<T>::_num_moved = 0;
template <class T>
int InstanceCounter<T>::_num_destroyed = 0;

#define EXPECT_INSTANCES(type, num_const, num_move, num_copy, num_inst)        \
    EXPECT_EQ((num_const), InstanceCounter<type>::getNumberOfConstructions()); \
    EXPECT_EQ((num_move), InstanceCounter<type>::getNumberOfMoves());          \
    EXPECT_EQ((num_copy), InstanceCounter<type>::getNumberOfCopies());         \
    EXPECT_EQ((num_inst), InstanceCounter<type>::getNumberOfInstances())

#define UPDATE_INSTANCES(type, num_const, num_move, num_copy, num_inst) \
    (num_const) = InstanceCounter<type>::getNumberOfConstructions();    \
    (num_move) = InstanceCounter<type>::getNumberOfMoves();             \
    (num_copy) = InstanceCounter<type>::getNumberOfCopies();            \
    (num_inst) = InstanceCounter<type>::getNumberOfInstances()
