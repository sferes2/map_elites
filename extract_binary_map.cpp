#include <iostream>
#include <algorithm>

#include "binary_map.hpp"

int main(int argc, char** argv)
{
    assert(argc == 2);
    binary_map::BinaryMap m = binary_map::load(argv[1]);
    std::vector<binary_map::Elem> v = m.elems;

    assert(v.size());
    std::copy(m.dims.begin(), m.dims.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    for (auto& x : v) {
        std::copy(x.pos.begin(), x.pos.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << x.fit << " ";
        std::copy(x.phen.begin(), x.phen.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;
    }
    return 0;
}
