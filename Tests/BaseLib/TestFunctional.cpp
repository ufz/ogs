/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include "BaseLib/Functional.h"
#include "Tests/InstanceCounter.h"

class A : private InstanceCounter<A>
{
public:
    A(const double value) : _value(value) {}
    A(A const& other) = default;
    A(A&& other) : InstanceCounter(std::move(other)), _value(other._value) {}
    A& operator=(A const& other) { _value = other._value; return *this; }
    A& operator=(A&& other) { _value = other._value; return *this; }

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

TEST(BaseLib, Functional)
{
    auto num_const = InstanceCounter<A>::getNumberOfConstructions();
    auto num_move = InstanceCounter<A>::getNumberOfMoves();
    auto num_copy = InstanceCounter<A>::getNumberOfCopies();
    auto num_inst = InstanceCounter<A>::getNumberOfInstances();

    // Base line: measure how many copies and moves
    // std::function<>(std::bind(...)) needs.
    A a_base(0.0);

    // move the object to std::bind()
    InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);
    std::function<void(A)> fct_mult(
        std::bind(&A::multiply, std::move(a_base), std::placeholders::_1));
    auto const num_copy_base_move =
        InstanceCounter<A>::getNumberOfCopies() - num_copy;
    auto const num_move_base_move =
        InstanceCounter<A>::getNumberOfMoves() - num_move;

    // call std::function using pass-by-value
    A a_base2(0.0);
    InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);
    fct_mult(a_base2);
    auto const num_copy_base_pass =
        InstanceCounter<A>::getNumberOfCopies() - num_copy;
    auto const num_move_base_pass =
        InstanceCounter<A>::getNumberOfMoves() - num_move;
    // end base line

    // self test
    InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);
    A a1(3.0);
    A a2(a1);
    EXPECT_INSTANCES(A, num_const+1, num_move, num_copy+1, num_inst+2);
    InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);

    auto f1_get = BaseLib::easyBind(&A::getValue, a1);
    EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
    EXPECT_EQ(3.0, f1_get());

    // check that really a reference is returned
    {
        auto f2_getRef = BaseLib::easyBind(&A::getValueRef, a2);
        auto& value_ref = f2_getRef();
        EXPECT_EQ(3.0, value_ref);
        value_ref = 4.0;
        EXPECT_EQ(4.0, a2.getValue());
    }
    EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);

    // test binding to pointers
    {
        A* ap = &a1;
        auto fp_get = BaseLib::easyBind(&A::getValue, ap);
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        EXPECT_EQ(3.0, fp_get());

        A const* apc = &a1;
        auto fpc_get = BaseLib::easyBind(&A::getValue, apc);
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        EXPECT_EQ(3.0, fpc_get());
    }
    EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);

    // check that referenced objects are not copied
    {
        A& a3 = a2;
        auto f3_get = BaseLib::easyBind(&A::getValue, a3);
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        EXPECT_EQ(4.0, f3_get());
    }
    {
        A const& a3 = a2;
        auto f3_get = BaseLib::easyBind(&A::getValue, a3);
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        EXPECT_EQ(4.0, f3_get());
    }
    EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);

    // temporaries must be moved
    {
        auto ftemp_get = BaseLib::easyBind(&A::getValue, A(5.0));

        EXPECT_INSTANCES(A, num_const + 1, num_move + num_move_base_move,
                         num_copy + num_copy_base_move, num_inst + 1);
        InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);

        EXPECT_EQ(5.0, ftemp_get());
    }
    // ftemp_get destroyed
    EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst-1);
    InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);

    // testing explicit move
    {
        A a_move(5.0);
        EXPECT_INSTANCES(A, num_const+1, num_move, num_copy, num_inst+1);
        InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);

        auto ftemp_get = BaseLib::easyBind(&A::getValue, std::move(a_move));

        EXPECT_INSTANCES(A, num_const, num_move + num_move_base_move,
                         num_copy + num_copy_base_move, num_inst+1);
        InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);

        EXPECT_EQ(5.0, ftemp_get());
    }
    // ftemp_get destroyed and a_move
    EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst-2);
    InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);

    // test binding a callable object
    {
        auto f1_op = BaseLib::easyBind(a1);
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        EXPECT_EQ(21.0, f1_op(7.0));
    }
    EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);

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
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        f1_add(a2);
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        EXPECT_EQ(7.0, f1_get());
    }

    // check that parameters passed by value are copied
    {
        auto f1_mult = BaseLib::easyBind(&A::multiply, a1);
        EXPECT_INSTANCES(A, num_const, num_move, num_copy, num_inst);
        f1_mult(a2);

        EXPECT_INSTANCES(A, num_const, num_move + num_move_base_pass,
                         num_copy + num_copy_base_pass, num_inst);
        InstanceCounter<A>::update(num_const, num_move, num_copy, num_inst);

        EXPECT_EQ(28.0, f1_get());
        EXPECT_EQ(4.0, a2.getValue());
    }
}
