#include <iostream>
#include <algorithm>

#include "binary_map.hpp"

int main(int argc, char **argv) {
  assert(argc == 2);
  std::vector<binary_map::Elem> v = binary_map::load(argv[1]);

  assert(v.size());
  for(auto &x : v) {
    std::copy(x.pos.begin(), x.pos.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << x.fit << " ";
    std::copy(x.phen.begin(), x.phen.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
  }
  return 0;
}
