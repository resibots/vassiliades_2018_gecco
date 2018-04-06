// Centroidal Voronoi Tesselation - Global Correlation
// Vassilis Vassiliades - December 2017

#ifndef CVT_GLOBAL_CORRELATION_HPP_
#define CVT_GLOBAL_CORRELATION_HPP_

#include <algorithm>
#include <limits>

#include <boost/array.hpp>
#include <boost/foreach.hpp>

#include <sferes/ea/ea.hpp>

#include <external/eigenmvn.h>

namespace sferes
{
namespace ea
{
// compute the matrix of distances
template <typename Phen>
struct _distance_f
{
  const Eigen::MatrixXf &_features;
  Eigen::MatrixXf &distances;

  ~_distance_f() {}
  _distance_f(const Eigen::MatrixXf &features, Eigen::MatrixXf &d)
      : _features(features), distances(d) {}

  _distance_f(const _distance_f &ev)
      : _features(ev._features), distances(ev.distances) {}

  void operator()(const parallel::range_t &r) const
  {
    for (size_t i = r.begin(); i != r.end(); ++i)
    {
      for (size_t j = 0; j < _features.rows(); ++j)
        distances(i, j) = (_features.row(i) - _features.row(j)).squaredNorm();
    }
  }
};

SFERES_EA(CVTGlobalCorrelation, Ea)
{
public:
  static constexpr size_t feature_dimensionality =
      Params::ea::feature_dimensionality;

  typedef boost::shared_ptr<Phen> indiv_t;
  typedef typename std::vector<indiv_t> pop_t;
  typedef boost::array<double, feature_dimensionality> point_t;

  static constexpr size_t number_of_clusters = Params::ea::number_of_clusters;
  std::vector<point_t> centroids;

  CVTGlobalCorrelation()
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

    // Convert the archive to an Eigen::MatrixXd
    std::pair<Eigen::MatrixXf, Eigen::MatrixXf> matrices = _get_eigen_matrices_from_pop(this->_pop);
    Eigen::MatrixXf &parameter_matrix = matrices.first;

    //Perform nearest-neighbor calculation in parameter space
    Eigen::MatrixXf &feature_matrix = matrices.first;

    // Compute pairwise distances and store them
    size_t archiveSize = this->_pop.size();

    Eigen::MatrixXf distances(archiveSize, archiveSize);
    _distance_f<Phen> f(feature_matrix, distances);
    parallel::init();
    parallel::p_for(parallel::range_t(0, archiveSize), f);

    // If we use the global correlation then there is no need to
    // take random points or nearest neighbors
    size_t k = archiveSize;

    // Fitting a Gaussian multivariate distribution in PARAMETER space
    Eigen::MatrixXf covar = Params::ea::alpha * _get_cov(parameter_matrix);

    for (size_t i = 0; i < Params::pop::size; ++i)
    {
      // Select a random parent (uniformly or according to some other probability distribution)
      int parent_index = _selection(this->_pop.size());

      // The mean is the PARAMETERS of the parent
      Eigen::VectorXf mean = parameter_matrix.row(parent_index);

      // Create a multivariate Gaussian distribution (Use the Cholesky decomposition)
      // Eigen::EigenMultivariateNormal<float> normX_cholesk(mean, covar);

      Eigen::VectorXf zero_vector = Eigen::VectorXf::Zero(Params::ea::dimensionality);
      Eigen::EigenMultivariateNormal<float> normX_cholesk(zero_vector, covar);

      // Sample 1 point from this distribution
      Eigen::MatrixXf samples = normX_cholesk.samples(1).transpose();

      // Add the mean
      samples.rowwise() += mean.adjoint();

      // Ensure that the sampled solutions are inside the bounds of the genotype space
      samples = _in_bounds(samples);

      // Get indiv_t individuals from Eigen
      pop_t sampled_individuals = _get_pop_from_eigen_matrix(samples);

      // Insert the individuals in the offspring array
      ptmp.insert(ptmp.end(), sampled_individuals.begin(), sampled_individuals.end());
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

  // We assume that bounds are in [0,1]
  Eigen::MatrixXf _in_bounds(const Eigen::MatrixXf &x)
  {
    Eigen::MatrixXf out = x;
    for (size_t i = 0; i < out.rows(); ++i)
    {
      for (size_t j = 0; j < out.cols(); ++j)
      {
        if (out(i, j) < 0.0)
          out(i, j) = 0.0;
        else if (out(i, j) > 1.0)
          out(i, j) = 1.0;
      }
    }
    return out;
  }

  // parameters to individuals
  pop_t _get_pop_from_eigen_matrix(const Eigen::MatrixXf &mat)
  {
    pop_t pop(mat.rows());

    for (size_t i = 0; i < pop.size(); ++i)
    {
      pop[i] = indiv_t(new Phen());

      // Set the data
      for (size_t j = 0; j < Params::ea::dimensionality; ++j)
        pop[i]->gen().data(j, mat(i, j));
    }

    return pop;
  }

  // individual parameters to matrix
  std::pair<Eigen::MatrixXf, Eigen::MatrixXf> _get_eigen_matrices_from_pop(const pop_t &pop)
  {
    assert(pop.size() > 0);

    Eigen::MatrixXf parameters = Eigen::MatrixXf::Zero(pop.size(), Params::ea::dimensionality);
    Eigen::MatrixXf features = Eigen::MatrixXf::Zero(pop.size(), feature_dimensionality);

    for (size_t i = 0; i < pop.size(); ++i)
    {
      for (size_t j = 0; j < parameters.cols(); ++j)
        parameters(i, j) = pop[i]->gen().data(j);

      for (size_t j = 0; j < features.cols(); ++j)
        features(i, j) = pop[i]->fit().desc(j);
    }

    return std::pair<Eigen::MatrixXf, Eigen::MatrixXf>(parameters, features);
  }

  // returns the covariance matrix
  // points are m x d
  Eigen::MatrixXf _get_cov(const Eigen::MatrixXf &points)
  {
    assert(points.rows() != 0 && points.cols() != 0);

    if (points.rows() == 1)
    {
      return Eigen::MatrixXf::Identity(points.cols(), points.cols());
    }
    else
    {
      Eigen::MatrixXf centered = points.rowwise() - points.colwise().mean();
      Eigen::MatrixXf cov = (centered.adjoint() * centered) / double(points.rows() - 1);

      return cov;
    }
  }

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

  int _selection(const size_t max_val)
  {
    return misc::rand<int>(0, max_val);
  }
};
}
}
#endif
