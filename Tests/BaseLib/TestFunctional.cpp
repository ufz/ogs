/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include "BaseLib/Functional.h"

class InstanceCounter
{
public:
    InstanceCounter() {
        ++_num_constructed;
    }
    InstanceCounter(InstanceCounter const&) {
        ++_num_copied;
    }
    InstanceCounter(InstanceCounter&&) {
        ++_num_moved;
    }
    ~InstanceCounter() {
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

private:
    static int _num_constructed;
    static int _num_copied;
    static int _num_moved;
    static int _num_destroyed;
};

int InstanceCounter::_num_constructed = 0;
int InstanceCounter::_num_copied = 0;
int InstanceCounter::_num_moved = 0;
int InstanceCounter::_num_destroyed = 0;

class A : public InstanceCounter
{
public:
    A(const double value) : _value(value) {}
    void add(A const& other) { _value += other._value; }
    // pass by value intended.
    void multiply(A other) {
        _value *= other._value;
        other._value = 0.0;
    }
    double getValue() const { return _value; }
    double& getValueRef() { return _value; }
    double operator()(double const x) { return _value*x; }
private:
    double _value;
};

#define EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst) \
    EXPECT_EQ((num_const), InstanceCounter::getNumberOfConstructions());    \
    EXPECT_EQ((num_move), InstanceCounter::getNumberOfMoves());             \
    EXPECT_EQ((num_copy), InstanceCounter::getNumberOfCopies());            \
    EXPECT_EQ((num_dest), InstanceCounter::getNumberOfDestructions());      \
    EXPECT_EQ((num_inst), InstanceCounter::getNumberOfInstances())

#define UPDATE_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst) \
    num_const = InstanceCounter::getNumberOfConstructions();              \
    num_move = InstanceCounter::getNumberOfMoves();                       \
    num_copy = InstanceCounter::getNumberOfCopies();                      \
    num_dest = InstanceCounter::getNumberOfDestructions();                \
    num_inst = InstanceCounter::getNumberOfInstances()

TEST(BaseLib, Functional)
{
    A a1(3.0);
    A a2(a1);

    auto num_const = InstanceCounter::getNumberOfConstructions();
    auto num_move = InstanceCounter::getNumberOfMoves();
    auto num_copy = InstanceCounter::getNumberOfCopies();
    auto num_dest = InstanceCounter::getNumberOfDestructions();
    auto num_inst = InstanceCounter::getNumberOfInstances();
    ASSERT_EQ(1, num_const);
    ASSERT_EQ(0, num_move);
    ASSERT_EQ(1, num_copy);
    ASSERT_EQ(0, num_dest);
    ASSERT_EQ(2, num_inst);

    auto f1_get = BaseLib::easyBind(&A::getValue, a1);
    EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
    EXPECT_EQ(3.0, f1_get());

    // check that really a reference is returned
    {
        auto f2_getRef = BaseLib::easyBind(&A::getValueRef, a2);
        auto& value_ref = f2_getRef();
        EXPECT_EQ(3.0, value_ref);
        value_ref = 4.0;
        EXPECT_EQ(4.0, a2.getValue());
    }
    EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);

    // test binding to pointers
    {
        A* ap = &a1;
        auto fp_get = BaseLib::easyBind(&A::getValue, ap);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        EXPECT_EQ(3.0, fp_get());

        A const* apc = &a1;
        auto fpc_get = BaseLib::easyBind(&A::getValue, apc);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        EXPECT_EQ(3.0, fpc_get());
    }
    EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);

    // check that referenced objects are not copied
    {
        A& a3 = a2;
        auto f3_get = BaseLib::easyBind(&A::getValue, a3);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        EXPECT_EQ(4.0, f3_get());
    }
    {
        A const& a3 = a2;
        auto f3_get = BaseLib::easyBind(&A::getValue, a3);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        EXPECT_EQ(4.0, f3_get());
    }
    EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);

    // temporaries must be moved
    {
        auto ftemp_get = BaseLib::easyBind(&A::getValue, A(5.0));

        EXPECT_EQ(num_const+1, InstanceCounter::getNumberOfConstructions());
        EXPECT_GE(num_move+2, InstanceCounter::getNumberOfMoves());
        EXPECT_EQ(num_copy, InstanceCounter::getNumberOfCopies());
        EXPECT_EQ(InstanceCounter::getNumberOfMoves(),
                  InstanceCounter::getNumberOfDestructions());
        EXPECT_EQ(num_inst+1, InstanceCounter::getNumberOfInstances());
        UPDATE_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);

        EXPECT_EQ(5.0, ftemp_get());
    }
    // ftemp_get destroyed
    EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest+1, num_inst-1);
    UPDATE_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);

    // test binding a callable object
    {
        auto f1_op = BaseLib::easyBind(a1);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        EXPECT_EQ(21.0, f1_op(7.0));
    }
    EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);

    // test binding a lambda
    {
        double value = 2.0;
        auto f_op = BaseLib::easyBind([&value](const double x) {
            value *= x;
            return value;
        });
        EXPECT_EQ(6.0, f_op(3.0));
        EXPECT_EQ(6.0, value);
    }

    // check that parameters passed by reference are not copied
    {
        auto f1_add = BaseLib::easyBind(&A::add, a1);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        f1_add(a2);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        EXPECT_EQ(7.0, f1_get());
    }

    // check that parameters passed by value are copied
    {
        auto f1_mult = BaseLib::easyBind(&A::multiply, a1);
        EXPECT_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);
        f1_mult(a2);

        EXPECT_EQ(num_const, InstanceCounter::getNumberOfConstructions());
        EXPECT_GE(num_move+2, InstanceCounter::getNumberOfMoves());
        EXPECT_EQ(num_copy+1, InstanceCounter::getNumberOfCopies());
        EXPECT_EQ(num_dest + InstanceCounter::getNumberOfMoves() +
                      InstanceCounter::getNumberOfCopies() - num_move -
                      num_copy,
                  InstanceCounter::getNumberOfDestructions());
        EXPECT_EQ(num_inst, InstanceCounter::getNumberOfInstances());
        UPDATE_INSTANCES(num_const, num_move, num_copy, num_dest, num_inst);

        EXPECT_EQ(28.0, f1_get());
        EXPECT_EQ(4.0, a2.getValue());
    }
}
