//[traverser1
#include <cassert>
#include <vector>
#include "mcrl2/core/traverser.h"

using namespace mcrl2;

struct int_traverser: public core::traverser<int_traverser>
{
  // The following lines are needed to make static polymorphism work.
  typedef core::traverser<int_traverser> super;
  using super::apply;
  using super::enter;
  using super::leave;

  int count = 0;

  // Handle integers.
  void apply(int i)
  {
    count += i;
  }
};

int main()
{
  std::vector<int> v = { 1, 4 };

  int_traverser f;
  f.apply(v);
  assert(f.count == 5);

  return 0;
}
//]