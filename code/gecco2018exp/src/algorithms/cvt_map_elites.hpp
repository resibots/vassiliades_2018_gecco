// Centroidal Voronoi Tesselation - MAP-Elites
// Vassilis Vassiliades - Inria, Nancy - 2018

#ifndef CVT_MAP_ELITES_HPP_
#define CVT_MAP_ELITES_HPP_

#include <algorithm>
#include <limits>

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/multi_array.hpp>

#include <sferes/ea/ea.hpp>
#include <sferes/fit/fitness.hpp>
#include <sferes/stc.hpp>

#include <Eigen/Core>

namespace sferes
{
namespace ea
{
SFERES_EA(CVTMapElites, Ea)
{
public:
  static constexpr size_t feature_dimensionality =
      Params::ea::feature_dimensionality;

  typedef boost::shared_ptr<Phen> indiv_t;
  typedef typename std::vector<indiv_t> pop_t;
  typedef boost::array<double, feature_dimensionality> point_t;

  static constexpr size_t number_of_clusters = Params::ea::number_of_clusters;
  std::vector<point_t> centroids;

  CVTMapElites()
  {
    _archive.resize(number_of_clusters);
    centroids = Params::ea::centroids;
  }

  void random_pop()
  {
    parallel::init();
    this->_pop.resize(Params::pop::init_size);

    BOOST_FOREACH (indiv_t &indiv, this->_pop)
    {
      indiv = indiv_t(new Phen());
      indiv->random();
    }

    this->_eval_pop(this->_pop, 0, this->_pop.size());

    std::vector<size_t> archive_indices = _get_archive_indices(this->_pop);

    for (size_t i = 0; i < this->_pop.size(); ++i)
      _add_to_archive(this->_pop[i], archive_indices[i]);
  }

  void epoch()
  {
    pop_t ptmp;

    this->_pop.clear();
    for (size_t i = 0; i < _archive.size(); ++i)
      if (_archive[i])
        this->_pop.push_back(_archive[i]);

    for (size_t i = 0; i < Params::pop::size; ++i)
    {
      size_t p1_index = _selection(this->_pop);
      indiv_t p1 = this->_pop[p1_index];
      size_t p2_index = p1_index;
      do{
        p2_index = _selection(this->_pop);
      } while(p1_index == p2_index);

      indiv_t p2 = this->_pop[p2_index];

      indiv_t i1, i2;
      p1->cross(p2, i1, i2);
      i1->mutate();
      ptmp.push_back(i1);
    }

    this->_eval_pop(ptmp, 0, ptmp.size());

    // Find the indices of the regions in which the descriptors map (nearest neighbor query)
    std::vector<size_t> offspring_centroid_indices = _get_archive_indices(ptmp);

    for (size_t i = 0; i < ptmp.size(); ++i)
      _add_to_archive(ptmp[i], offspring_centroid_indices[i]);
  }

  const pop_t &archive() const { return _archive; }

  template <typename I>
  point_t get_point(const I &indiv) const
  {
    return _get_point(indiv);
  }

protected:
  pop_t _archive;

  std::vector<size_t> _get_archive_indices(pop_t & pop)
  {
    size_t pop_size = pop.size();
    std::vector<size_t> archive_indices(pop_size, -1);

    tbb::parallel_for(size_t(0), pop_size, size_t(1), [&](size_t i) {
      indiv_t &indiv = pop[i];
      point_t p = _get_point(indiv);
      double min_dist = std::numeric_limits<double>::max();
      size_t archive_index = -1;

      // Find the closest centroid
      for (size_t i = 0; i < centroids.size(); ++i)
      {
        double dist = _calc_dist(centroids[i], p);

        if (dist < min_dist)
        {
          min_dist = dist;
          archive_index = i;
        }

        // Since the minimum distance cannot be less than 0
        // we could accelerate computation by breaking
        if (min_dist == 0.0)
          break;
      }

      // Store the archive index of the i'th individual
      archive_indices[i] = archive_index;
    });

    return archive_indices;
  }

  bool _add_to_archive(indiv_t i1, const size_t archive_index)
  {
    if (i1->fit().dead())
      return false;

    // If the archive is empty or the stored individual is less fit
    if (!_archive[archive_index] ||
        i1->fit().value() > _archive[archive_index]->fit().value())
    {
      _archive[archive_index] = i1;
      return true;
    }

    return false;
  }

  // Euclidean distance
  double _calc_dist(const point_t &p1, const point_t &p2)
  {
    double dist = 0.0;

    for (size_t i = 0; i < feature_dimensionality; ++i)
      dist += pow(p1[i] - p2[i], 2);

    return sqrt(dist);
  }

  template <typename I>
  point_t _get_point(const I &indiv) const
  {
    point_t p;
    for (size_t i = 0; i < feature_dimensionality; ++i)
      p[i] = indiv->fit().desc()[i];

    return p;
  }

  size_t _selection(const pop_t &pop)
  {
    return misc::rand<int>(0, pop.size());
  }
};
}
}
#endif
