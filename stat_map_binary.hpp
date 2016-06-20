#ifndef STAT_MAP_BINARY_HPP_
#define STAT_MAP_BINARY_HPP_

#include <numeric>
#include <boost/multi_array.hpp>
#include <sferes/stat/stat.hpp>

#include "binary_map.hpp"

namespace sferes {
    namespace stat {

        // - write in generic binary format
        // - work only for evofloat & sampled !
        // (the binary map in sferes form is in gen_* if you have StatMap)
        SFERES_STAT(MapBinary, Stat){
            public :
                template <typename E>
                void refresh(const E& ea){
                    if (ea.gen() % Params::pop::dump_period != 0) return;
        std::string fname = ea.res_dir() + "/"
            + "map_" + boost::lexical_cast<std::string>(ea.gen())
            + std::string(".bin");
        auto archive = _multi_array_to_vector(ea.archive(), ea);
        binary_map::write(archive, fname);
    }

protected:
    // convert a (possibly sparse) boost::multi_array of phenotypes
    // to a list of vectors that can be used in other software
    template <typename A, typename EA>
    binary_map::BinaryMap _multi_array_to_vector(const A& array, const EA& ea)
    {
        binary_map::BinaryMap map;
        map.dims.clear();
        map.elems.clear();
        for (size_t dim = 0; dim < array.num_dimensions(); ++dim)
            map.dims.push_back(array.shape()[dim]);
        size_t offset = 0;
        for (auto i = array.data(); i < (array.data() + array.num_elements()); ++i) {
            if (*i) {
                auto posinarray = ea.getindexarray(array, i);
                std::vector<int> pos;
                for (size_t dim = 0; dim < array.num_dimensions(); ++dim) {
                    pos.push_back(posinarray[dim]);
                }

                auto phen = array(posinarray);
                float fit = phen->fit().value();
                std::vector<float> p;
                for (size_t j = 0; j < phen->size(); ++j)
                    p.push_back(phen->data(j));
                map.elems.push_back(binary_map::Elem{pos, p, fit});
            }
            ++offset;
        }
        return map;
    }
};
}
}

#endif
