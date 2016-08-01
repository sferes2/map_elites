#ifndef STAT_MAP_BINARY_HPP_
#define STAT_MAP_BINARY_HPP_

#include <boost/multi_array.hpp>
#include <numeric>
#include <sferes/stat/stat.hpp>

#include "binary_map.hpp"

namespace sferes {
namespace stat {

// - write in generic binary format
// - work only for evofloat & sampled !
// (the binary map in sferes form is in gen_* if you have StatMap)
SFERES_STAT(MapBinary,
            Stat){public : template <typename E> void refresh(const E &ea){
                      if (ea.gen() % Params::pop::dump_period != 0) return;
std::string fname = ea.res_dir() + "/" + "map_" +
                    boost::lexical_cast<std::string>(ea.gen()) +
                    std::string(".bin");
auto archive = _multi_array_to_vector(ea.archive(), ea);
binary_map::write(archive, fname);
}

protected:
/*      template<typename A>
  long int get_index(const A& m, const typename A::element* e, size_t dir) {
  int offset = e - m.origin();
  return (offset / m.strides()[dir] % m.shape()[dir] +  m.index_bases()[dir]);
}

template<typename A>
  std::array<typename A::index, A::dimensionality> get_indexarray(const A& m,
const typename A::element* e) {
  std::array<typename A::index, A::dimensionality> _index;
  for (size_t dir = 0; dir < m.num_dimensions(); dir++) {
    _index[dir] = get_index(m, e, dir);
  }
  return _index;
}
*/

// convert a (possibly sparse) boost::multi_array of phenotypes
// to a list of vectors that can be used in other software
template <typename A, typename EA>
std::vector<binary_map::Elem> _multi_array_to_vector(const A &array,
                                                     const EA &ea) {
  std::vector<binary_map::Elem> res;
  size_t offset = 0;
  for (auto i = array.data(); i < (array.data() + array.num_elements()); ++i) {
    if (*i) {
      auto posinarray = ea.getindexarray(array, i);
      std::vector<int> pos; // double?
      for (size_t dim = 0; dim < array.num_dimensions(); ++dim)
        pos.push_back(posinarray[dim]); // / (double) m.shape()[dim]);

      auto phen = array(posinarray);
      float fit = phen->fit().value();
      std::vector<float> p;
      for (size_t j = 0; j < phen->size(); ++j)
        p.push_back(phen->data(j));
      res.push_back(binary_map::Elem{pos, p, fit});
    }
    ++offset;
  }
  return res;
}
};
}
}

#endif
