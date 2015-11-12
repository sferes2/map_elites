#ifndef STAT_MAP_HPP_
#define STAT_MAP_HPP_

#include <numeric>
#include <boost/multi_array.hpp>
#include <sferes/stat/stat.hpp>

#define MAP_WRITE_PARENTS

namespace sferes {
  namespace stat {
    SFERES_STAT(Map, Stat) {
    public:
      typedef boost::shared_ptr<Phen> phen_t;
      typedef boost::multi_array<phen_t, Params::ea::behav_dim> array_t;
      typedef boost::array<float, Params::ea::behav_dim> point_t;
      typedef boost::array<typename array_t::index, Params::ea::behav_dim> behav_index_t;

      size_t behav_dim;
      behav_index_t behav_shape;
      behav_index_t behav_strides;
      behav_index_t behav_indexbase;

      Map() : behav_dim(Params::ea::behav_dim) {
        for(size_t i = 0; i < Params::ea::behav_shape_size(); ++i)
          behav_shape[i] = Params::ea::behav_shape(i);
      }

      template<typename E>
      void refresh(const E& ea) {
        _archive.clear();

        for(size_t i = 0; i < behav_dim; ++i) {
          assert(ea.archive().shape()[i] == behav_shape[i]);
          behav_strides[i] = ea.archive().strides()[i];
          behav_indexbase[i] = ea.archive().index_bases()[i];
        }

        for(const phen_t* i = ea.archive().data(); i < (ea.archive().data() + ea.archive().num_elements()); ++i) {
          phen_t p = *i;
          _archive.push_back(p);
        }

        this->_create_log_file(ea, "progress_archive.dat");
        _write_progress(ea, *this->_log_file);

        if (ea.gen() % Params::pop::dump_period == 0) {
          _write_archive(ea.archive(), std::string("archive_"), ea);
#ifdef MAP_WRITE_PARENTS
          _write_parents(ea.archive(), ea.parents(), std::string("parents_"), ea);
#endif
        }
      }


      void show(std::ostream& os, size_t k) {
        std::cerr << "loading " ;
        for(size_t i = 0; i < behav_dim; ++i)
          std::cerr << (k / behav_strides[i] % behav_shape[i] +  behav_indexbase[i]) << ",";
        std::cerr << std::endl;


        if (_archive[k]) {
          _archive[k]->develop();
          _archive[k]->show(os);
          _archive[k]->fit().set_mode(fit::mode::view);
          _archive[k]->fit().eval(*_archive[k]);
        } else
          std::cerr << "Warning, no point here" << std::endl;
      }

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_NVP(_archive);
        ar & BOOST_SERIALIZATION_NVP(behav_dim);
        ar & BOOST_SERIALIZATION_NVP(behav_shape);
        ar & BOOST_SERIALIZATION_NVP(behav_strides);
        ar & BOOST_SERIALIZATION_NVP(behav_indexbase);
      }



    protected:
      std::vector<phen_t> _archive;
      //int _xs, _ys;

      template<typename EA>
      void _write_parents(const array_t& array,
                          const array_t& p_array,
                          const std::string& prefix,
                          const EA& ea) const {
        std::cout << "writing..." << prefix << ea.gen() << std::endl;
        std::string fname =  ea.res_dir() + "/"
                             + prefix
                             + boost::lexical_cast<
                             std::string>(ea.gen())
                             + std::string(".dat");
        std::ofstream ofs(fname.c_str());

        for(const phen_t* i = array.data(); i < (array.data() + array.num_elements()); ++i) {
          if (*i) {
            behav_index_t idx = ea.getindexarray(array, i);
            assert(array(idx)->fit().value() == (*i)->fit().value());
            if(p_array(idx)) {
              for(size_t dim = 0; dim < behav_dim; ++dim)
                ofs << idx[dim] / (float) behav_shape[dim] << " ";
              ofs << " " << p_array(idx)->fit().value() << " " ;

              point_t p = ea.get_point(p_array(idx));
              behav_index_t posinparent;
              for(size_t dim = 0; dim < behav_dim; ++dim) {
                posinparent[dim] = round(p[dim] * behav_shape[dim]);
                ofs << posinparent[dim] / (float) behav_shape[dim] << " ";
              }
              ofs << " " << array(idx)->fit().value() << std::endl;
            }
          }
        }
      }

      template<typename EA>
      void _write_archive(const array_t& array,
                          const std::string& prefix,
                          const EA& ea) const {
        std::cout << "writing..." << prefix << ea.gen() << std::endl;
        std::string fname = ea.res_dir() + "/"
                            + prefix
                            + boost::lexical_cast<
                            std::string>(ea.gen())
                            + std::string(".dat");

        std::ofstream ofs(fname.c_str());

        size_t offset = 0;
        for(const phen_t* i = array.data(); i < (array.data() + array.num_elements()); ++i) {
          if (*i) {
            behav_index_t posinarray = ea.getindexarray(array, i);
            assert(array(posinarray)->fit().value() == (*i)->fit().value());

            ofs << offset << "    ";
            for(size_t dim = 0; dim < behav_dim; ++dim)
              ofs << posinarray[dim] / (float) behav_shape[dim] << " ";
            ofs << " " << array(posinarray)->fit().value() << std::endl;
          }
          ++offset;
        }

      }

      template<typename EA>
      void _write_progress(const EA& ea, std::ofstream& ofs) const {

        const array_t& array = ea.archive();

        size_t archive_size = 0;
        float archive_mean = 0.0f;
        float archive_max = -1e10;
        float mean_dist_center = 0.0f;

        for(const phen_t* i = array.data(); i < (array.data() + array.num_elements()); ++i) {
          if(*i) {
            phen_t p = *i;
            archive_size++;
            archive_mean += p->fit().value();

            if(archive_max < p->fit().value())
              archive_max = p->fit().value();

            float dist = 0.0f;
            for(size_t i = 0; i < Params::ea::behav_shape_size(); ++i) {
              assert(p->fit().desc()[i] >= 0.0f && p->fit().desc()[i] <= 1.0f);

              float diff = p->fit().desc()[i] -
                           (float)round(p->fit().desc()[i] * (float)(Params::ea::behav_shape(i)-1)) / (float)(Params::ea::behav_shape(i) - 1);

              dist += diff * diff;
            }
            mean_dist_center+=sqrtf(dist);
          }
        }

        archive_mean /= archive_size;
        mean_dist_center /= archive_size;
        ofs << ea.gen() << " " << archive_size << " " << archive_mean << " " << archive_max << " " << mean_dist_center << std::endl;
      }


    };
  }
}

#endif
