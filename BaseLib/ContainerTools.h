/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <type_traits>

namespace BaseLib
{
namespace detail
{

template <typename Element>
struct PolymorphicRandomAccessContainerViewInterface;

template <typename Element>
struct PolymorphicRandomAccessContainerViewIterator
{
    using value_type = Element;
    using difference_type = std::ptrdiff_t;

    PolymorphicRandomAccessContainerViewIterator() = default;

    PolymorphicRandomAccessContainerViewIterator(
        std::size_t n,
        PolymorphicRandomAccessContainerViewInterface<Element> const* const
            view) noexcept
        : n_{n}, view_{view}
    {
    }

    PolymorphicRandomAccessContainerViewIterator& operator++() noexcept
    {
        ++n_;
        return *this;
    }
    PolymorphicRandomAccessContainerViewIterator operator++(int) noexcept
    {
        auto copy{*this};
        operator++();
        return copy;
    }

    PolymorphicRandomAccessContainerViewIterator& operator--() noexcept
    {
        --n_;
        return *this;
    }
    PolymorphicRandomAccessContainerViewIterator operator--(int) noexcept
    {
        auto copy{*this};
        operator--();
        return copy;
    }

    PolymorphicRandomAccessContainerViewIterator& operator+=(
        difference_type const increment) noexcept
    {
        n_ += increment;
        return *this;
    }
    PolymorphicRandomAccessContainerViewIterator& operator-=(
        difference_type const decrement) noexcept
    {
        n_ -= decrement;
        return *this;
    }

    [[nodiscard]] friend PolymorphicRandomAccessContainerViewIterator operator+(
        PolymorphicRandomAccessContainerViewIterator iterator,
        difference_type const increment) noexcept
    {
        iterator.n_ += increment;
        return iterator;
    }

    [[nodiscard]] friend PolymorphicRandomAccessContainerViewIterator operator+(
        difference_type const increment,
        PolymorphicRandomAccessContainerViewIterator iterator) noexcept
    {
        iterator.n_ += increment;
        return iterator;
    }

    [[nodiscard]] difference_type operator-(
        PolymorphicRandomAccessContainerViewIterator const other) const noexcept
    {
        assert(view_ == other.view_);
        return n_ - other.n_;
    }

    [[nodiscard]] friend PolymorphicRandomAccessContainerViewIterator operator-(
        PolymorphicRandomAccessContainerViewIterator iterator,
        difference_type const decrement) noexcept
    {
        iterator.n_ -= decrement;
        return iterator;
    }

    [[nodiscard]] Element& operator*() const { return (*view_)[n_]; }
    [[nodiscard]] Element& operator[](difference_type n) const
    {
        return (*view_)[n_ + n];
    }

    auto operator<=>(
        PolymorphicRandomAccessContainerViewIterator const&) const = default;

private:
    std::size_t n_ = 0;
    PolymorphicRandomAccessContainerViewInterface<Element> const* view_ =
        nullptr;
};

template <typename Element>
struct PolymorphicRandomAccessContainerViewInterface
{
    [[nodiscard]] PolymorphicRandomAccessContainerViewIterator<Element> begin()
        const noexcept
    {
        return {0, this};
    }
    [[nodiscard]] PolymorphicRandomAccessContainerViewIterator<Element> end()
        const noexcept
    {
        return {size(), this};
    }

    [[nodiscard]] virtual std::size_t size() const noexcept = 0;

    [[nodiscard]] virtual Element& operator[](std::size_t) const = 0;

    virtual ~PolymorphicRandomAccessContainerViewInterface() = default;
};

template <typename Element, typename Container>
struct CovariantRandomAccessContainerView final
    : PolymorphicRandomAccessContainerViewInterface<Element>
{
    explicit CovariantRandomAccessContainerView(Container& container)
        : container_{container}
    {
    }

    Element& operator[](std::size_t index) const override
    {
        if constexpr (requires {
                          typename Container::value_type::element_type;
                      })
        {
            static_assert(
                std::derived_from<typename Container::value_type::element_type,
                                  Element>);

            return *container_[index];
        }
        else
        {
            static_assert(
                std::derived_from<typename Container::value_type, Element>);

            // this restriction is necessary because of the constness
            // propagation of std::vector.
            static_assert(
                (!std::is_const_v<Container>) || std::is_const_v<Element>,
                "If the element type is non-const, the container must be "
                "non-const, too.");

            return container_[index];
        }
    }

    std::size_t size() const noexcept override { return container_.size(); }

private:
    Container& container_;
};
}  // namespace detail

/**
 * Provides access to the (possibly upcasted) elements of a random access
 * container.
 *
 * Example:
 *
 * \code{.cpp}
 * std::vector<std::unique_ptr<Pet>> pets;
 * pets.push_back(std::make_unique<Budgie>());
 * pets.push_back(std::make_unique<Dog>());
 *
 * // Pet is derived from Animal
 * auto f = [](BaseLib::PolymorphicRandomAccessContainerView<Animal> const&
 *              animals_)
 * {
 *  ASSERT_EQ(2, animals_.size());
 *
 *  for (auto& animal : animals_)
 *  {
 *      animal.color = "red";
 *  }
 * };
 *
 * f(pets);  // pass the std::vector<std::unique_ptr<Pet>> to the function
 *           // taking the view
 * \endcode
 *
 * For further usage examples take a look into the unit tests.
 */
template <typename Element>
struct PolymorphicRandomAccessContainerView
{
    PolymorphicRandomAccessContainerView(
        PolymorphicRandomAccessContainerView const&) = delete;
    PolymorphicRandomAccessContainerView(
        PolymorphicRandomAccessContainerView&&) = default;

    template <typename Container>
    explicit(false)
        // cppcheck-suppress noExplicitConstructor
        PolymorphicRandomAccessContainerView(Container& container)
        : impl_{std::make_unique<
              detail::CovariantRandomAccessContainerView<Element, Container>>(
              container)}
    {
    }

    template <typename Container>
    explicit(false)
        // cppcheck-suppress noExplicitConstructor
        PolymorphicRandomAccessContainerView(Container const& container)
        : impl_{std::make_unique<
              detail::CovariantRandomAccessContainerView<Element,
                                                         Container const>>(
              container)}
    {
    }

    PolymorphicRandomAccessContainerView& operator=(
        PolymorphicRandomAccessContainerView const&) = delete;
    PolymorphicRandomAccessContainerView& operator=(
        PolymorphicRandomAccessContainerView&&) = default;

    [[nodiscard]] detail::PolymorphicRandomAccessContainerViewIterator<Element>
    begin() const noexcept
    {
        return impl_->begin();
    }
    [[nodiscard]] detail::PolymorphicRandomAccessContainerViewIterator<Element>
    end() const noexcept
    {
        return impl_->end();
    }

    [[nodiscard]] std::size_t size() const noexcept { return impl_->size(); }

    [[nodiscard]] Element& operator[](std::size_t n) const
    {
        return (*impl_)[n];
    }

private:
    std::unique_ptr<
        detail::PolymorphicRandomAccessContainerViewInterface<Element>>
        impl_;
};
}  // namespace BaseLib
