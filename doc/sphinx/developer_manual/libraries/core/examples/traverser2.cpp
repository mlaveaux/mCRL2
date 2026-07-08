//[traverser2
#include <cassert>
#include <iostream>
#include <vector>
#include "mcrl2/core/traverser.h"

using namespace mcrl2;

struct A
{
  int i;
  std::vector<int> v;
  
  A(int i0, const std::vector<int>& v0)
    : i(i0), v(v0)
  {}
};

template <typename Derived>
struct int_traverser: public core::traverser<Derived>
{
  typedef core::traverser<Derived> super;
  using super::apply;
  using super::enter;
  using super::leave;

  // Handle integers.
  void apply(int i)
  {
    std::cout << "i = " << i << '\n';
  }
};

template <typename Derived>
struct A_traverser: public int_traverser<Derived>
{
  typedef int_traverser<Derived> super;
  using super::apply;
  using super::enter;
  using super::leave;

  // Handle A.
  void apply(const A& a)
  {
    // The static_cast<Derived&>(*this) makes it possible to override functions.
    static_cast<Derived&>(*this).apply(a.i);
    static_cast<Derived&>(*this).apply(a.v);
  }
};

template <typename Derived>
struct B_traverser: public A_traverser<Derived>
{
  typedef A_traverser<Derived> super;
  using super::apply;
  using super::enter;
  using super::leave;

  // Override the handler for integers.
  void apply(int i)
  {
    std::cout << "2*i = " << 2*i << '\n';
  }
};

int main()
{
  int i = 3;
  std::vector<int> v;
  v.push_back(1);
  v.push_back(4);
  A a(i, v);

  // Apply A_traverser to a. This gives the following output:
  //
  // i = 3
  // i = 1
  // i = 4
  //
  core::apply_traverser<A_traverser>().apply(a);

  // Apply B_traverser to a. This gives the following output:
  //
  // 2*i = 6
  // 2*i = 2
  // 2*i = 8
  //
  core::apply_traverser<B_traverser>().apply(a);

  return 0;
}
//]