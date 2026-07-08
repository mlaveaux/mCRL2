//[builder1
#include <cassert>
#include <vector>
#include "mcrl2/core/builder.h"

using namespace mcrl2;

struct int_builder: public core::builder<int_builder>
{
  // The following lines are needed to make static polymorphism work.
  typedef core::builder<int_builder> super;
  using super::enter;
  using super::leave;
  using super::update;

  // Handle integers.
  void apply(int& result, int i)
  {
    result = i + 1;
  }
};

int main()
{
  std::vector<int> v = { 1, 4 };

  int_builder f;
  f.update(v);
  assert(v[0] == 2);
  assert(v[1] == 5);

  return 0;
}
//]