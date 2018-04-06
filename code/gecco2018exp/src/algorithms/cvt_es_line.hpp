// Centroidal Voronoi Tesselation - Evolution Strategies + Line
// Vassilis Vassiliades - December 2017

#ifndef CVT_ES_LINE_HPP_
#define CVT_ES_LINE_HPP_

#include <Eigen/Core>
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

#include <random>

#include "archive_elem_es_line.hpp"

namespace sferes
{
namespace ea
{
SFERES_EA(CVTESLine, Ea)
{
  public:
    typedef boost::shared_ptr<archive_elem_es_line<Params, Phen>> archive_elem_ptr_t;
    typedef std::vector<archive_elem_ptr_t> archive_t;

    static constexpr size_t feature_dimensionality =
        Params::ea::feature_dimensionality;

    static constexpr size_t dimensionality =
        Params::ea::dimensionality;

    typedef boost::shared_ptr<Phen> indiv_t;
    typedef typename std::vector<indiv_t> pop_t;
    typedef boost::array<double, feature_dimensionality> point_t;

    static constexpr size_t number_of_clusters = Params::ea::number_of_clusters;
    std::vector<point_t> centroids;

    CVTESLine()
    {
        _archive.resize(number_of_clusters);
        centroids = Params::ea::centroids;

        _mersenne_engine = std::mt19937(_rd()); // Standard mersenne_twister_engine seeded with rd()
        _normal_dis = std::normal_distribution<>(0, 1);
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
        {
            archive_elem_ptr_t elem = archive_elem_ptr_t(new archive_elem_es_line<Params, Phen>(this->_pop[i]));
            _add_to_archive(elem, archive_indices[i]);
        }

        archive_t _temp_archive;
        for (size_t i = 0; i < _archive.size(); ++i)
            if (_archive[i])
                _temp_archive.push_back(_archive[i]);
    }

    void epoch()
    {
        archive_t _temp_archive;
        pop_t ptmp;

        // Get the elements
        for (size_t i = 0; i < _archive.size(); ++i)
            if (_archive[i])
                _temp_archive.push_back(_archive[i]);

        archive_t offspring_archive;

        for (size_t i = 0; i < Params::pop::size; ++i)
        {
            // Select a parent from a (filled) region
            size_t p1_index = _selection(_temp_archive);
            archive_elem_ptr_t parent1 = _temp_archive[p1_index];

            // Deep copy of the parent (indiv and strategy parameters)
            archive_elem_ptr_t offspring = archive_elem_ptr_t(new archive_elem_es_line<Params, Phen>(*parent1));

            // Get an Eigen object of the offspring genotype (to manipulate it easier)
            Eigen::VectorXd sol = _get_genotype(offspring);

            // Sample from normal distribution
            Eigen::VectorXd z_sigma_isotropic = _sample_gaussian(1);
            Eigen::VectorXd z_sigma_directed = _sample_gaussian(1);

            // Adapt the sigma
            offspring->sigma_isotropic = offspring->sigma_isotropic * exp(Params::ea::tau_isotropic * z_sigma_isotropic(0));
            offspring->sigma_directed = offspring->sigma_directed * exp(Params::ea::tau_directed * z_sigma_directed(0));

            // Create the new solution
            Eigen::VectorXd z_isotropic = _sample_gaussian(dimensionality);
            Eigen::VectorXd z_directed = _sample_gaussian(1);

            // Find the direction
            size_t p2_index = p1_index;
            archive_elem_ptr_t parent2;

            // don't allow selecting the same parent
            do {
                p2_index = _selection(_temp_archive);
                parent2 = _temp_archive[p2_index];
            } while(p1_index == p2_index);

            Eigen::VectorXd gen2 = _get_genotype(parent2);
            Eigen::VectorXd direction = gen2 - sol;
            double direction_norm = direction.norm();

            // Handle the scaling by distance
            if (Params::ea::isotropic_scaled_by_distance)
                z_isotropic = z_isotropic / direction_norm;

            if (!Params::ea::directed_scaled_by_distance)
                direction = direction / direction_norm;

            Eigen::VectorXd new_sol = sol + offspring->sigma_isotropic * z_isotropic + offspring->sigma_directed * direction * z_directed;

            // Ensure that the new solution is in the bounds of the genotype space
            new_sol = _in_bounds(new_sol);

            // Get back the genotype from the Eigen object
            for (size_t j = 0; j < new_sol.size(); ++j)
                offspring->indiv->gen().data(j, (float)new_sol(j));

            // Store the individual to be evaluated
            ptmp.push_back(offspring->indiv);

            // Store the offspring to add it to the archive
            offspring_archive.push_back(offspring);
        }

        // Evaluation
        this->_eval_pop(ptmp, 0, ptmp.size());

        // Find the indices of the regions in which the descriptors map (nearest neighbor query)
        std::vector<size_t> offspring_centroid_indices = _get_archive_indices(ptmp);

        for (size_t i = 0; i < offspring_archive.size(); ++i)
            _add_to_archive(offspring_archive[i], offspring_centroid_indices[i]);
    }

    const std::vector<archive_elem_ptr_t> &archive() const { return _archive; }

    template <typename I>
    point_t get_point(const I &indiv) const
    {
        return _get_point(indiv);
    }

  protected:
    archive_t _archive;

    std::random_device _rd;        // Will be used to obtain a seed for the random number engine
    std::mt19937 _mersenne_engine; // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> _normal_dis;

    Eigen::VectorXd _get_genotype(const archive_elem_ptr_t &elem)
    {
        Eigen::VectorXd genotype = Eigen::VectorXd::Zero(dimensionality);
        for (size_t j = 0; j < dimensionality; ++j)
            genotype(j) = elem->indiv->gen().data(j);
        return genotype;
    }

    Eigen::VectorXd _sample_gaussian(const size_t dim)
    {
        std::vector<double> x(dim, 0.0);
        std::generate(x.begin(), x.end(), [this]() { return this->_normal_dis(_mersenne_engine); });

        Eigen::VectorXd out(dim);
        out = Eigen::VectorXd::Map(x.data(), x.size());

        return out;
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

    bool _add_to_archive(archive_elem_ptr_t elem, const size_t archive_index)
    {
        if (elem->indiv->fit().dead())
            return false;

        if (!_archive[archive_index])
        {
            _archive[archive_index] = elem;
            return true;
        }
        else if (elem->indiv->fit().value() > _archive[archive_index]->indiv->fit().value())
        {
            elem->parent_index = archive_index;
            _archive[archive_index] = elem;
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

    // archive_elem_ptr_t _selection(const archive_t &a)
    // {
    //     int index = misc::rand<int>(0, a.size());
    //     return a[index];
    // }

    size_t _selection(const archive_t &a)
    {
        return misc::rand<int>(0, a.size());
    }

    // We assume that bounds are in [0,1]
    Eigen::VectorXd _in_bounds(const Eigen::VectorXd &x)
    {
        Eigen::VectorXd out = x;
        for (size_t i = 0; i < out.size(); ++i)
        {
            if (out(i) < 0.0)
                out(i) = 0.0;
            else if (out(i) > 1.0)
                out(i) = 1.0;
        }
        return out;
    }
};
}
}
#endif
