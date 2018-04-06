// Vassilis Vassiliades - Inria, Nancy - 2018

#include <iostream>
#include <sferes/gen/evo_float.hpp>
#include <sferes/modif/dummy.hpp>
#include <sferes/phen/parameters.hpp>
#include <sferes/run.hpp>

#define NO_MPI
#include <sferes/eval/parallel.hpp>

#if defined(VARIATIONGC)
#include "algorithms/cvt_global_correlation.hpp"
#elif defined(VARIATIONSBX)
#include "algorithms/cvt_map_elites.hpp"
#else
#include "algorithms/cvt_es_line.hpp"
#endif

#include "algorithms/fit_map.hpp"
#include "algorithms/stat_map.hpp"
#include "algorithms/archive_fit.hpp"

using namespace sferes;
using namespace sferes::gen::evo_float;

struct Params
{
    struct ea
    {
        SFERES_CONST size_t number_of_clusters = 10000;

        SFERES_CONST size_t dimensionality = 100;
        SFERES_CONST size_t feature_dimensionality = 2;
        typedef boost::array<double, feature_dimensionality> point_t;
        static std::vector<point_t> centroids;

        // No self-adaptation for the directed part
        SFERES_CONST double tau_directed = 0.0;

#ifdef VARIATIONGC
        SFERES_CONST double alpha = 0.1;
#elif VARIATIONISO
        SFERES_CONST double sigma_isotropic = 0.1;
        SFERES_CONST double sigma_directed = 0.0;

        // No self-adaptation
        SFERES_CONST double tau_isotropic = 0.0;

        // No scaling by distance
        SFERES_CONST bool isotropic_scaled_by_distance = false;
        SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONISOSA
        SFERES_CONST double sigma_isotropic = 0.1; // this acts as an initial value
        SFERES_CONST double sigma_directed = 0.0;

        // Self-adaptation
        SFERES_CONST double tau_isotropic = 1.0 / sqrt(2.0 * dimensionality);

        // No scaling by distance
        SFERES_CONST bool isotropic_scaled_by_distance = false;
        SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONISODD
        SFERES_CONST double sigma_isotropic = 0.05;
        SFERES_CONST double sigma_directed = 0.0;

        // No self-adaptation
        SFERES_CONST double tau_isotropic = 0.0;

        // Scaling by distance for the isotropic
        SFERES_CONST bool isotropic_scaled_by_distance = true;
        SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONISOLINEDD
        SFERES_CONST double sigma_isotropic = 0.01;
        SFERES_CONST double sigma_directed = 0.2;

        // No self-adaptation
        SFERES_CONST double tau_isotropic = 0.0;

        // Scaling by distance for the directed
        SFERES_CONST bool isotropic_scaled_by_distance = false;
        SFERES_CONST bool directed_scaled_by_distance = true;
#elif VARIATIONLINE
    SFERES_CONST double sigma_isotropic = 0.0;
    SFERES_CONST double sigma_directed = 0.2;

    // No self-adaptation
    SFERES_CONST double tau_isotropic = 0.0;

    // No scaling by distance
    SFERES_CONST bool isotropic_scaled_by_distance = false;
    SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONLINEDD
    SFERES_CONST double sigma_isotropic = 0.0;
    SFERES_CONST double sigma_directed = 0.2;

    // No self-adaptation
    SFERES_CONST double tau_isotropic = 0.0;

    // Scaling by distance for the directed
    SFERES_CONST bool isotropic_scaled_by_distance = false;
    SFERES_CONST bool directed_scaled_by_distance = true;
#endif
    };

    struct pop
    {
        SFERES_CONST size_t init_size = 100;
        SFERES_CONST size_t size = 100;
        SFERES_CONST size_t nb_gen = 999;
        SFERES_CONST size_t dump_period = 500;
    };

    struct evo_float
    {
        SFERES_CONST float cross_rate = 1.0f;

#if defined(VARIATIONSBX)
        SFERES_CONST float eta_c = 10.0f;
        SFERES_CONST cross_over_t cross_over_type = sbx;
#else
        SFERES_CONST cross_over_t cross_over_type = no_cross_over;
#endif
        SFERES_CONST float mutation_rate = 0.0f;
        SFERES_CONST float sigma = 0.0f;
        SFERES_CONST mutation_t mutation_type = gaussian;
    };

    struct parameters
    {
        SFERES_CONST float min = -5.0;
        SFERES_CONST float max = 5.0;
    };
};

typedef Params::ea::point_t point_t;

// define the centroids
std::vector<point_t> Params::ea::centroids;

template <typename Point>
std::vector<Point> load_centroids(const std::string &centroids_filename)
{
    std::vector<Point> centroids;

    std::ifstream fin(centroids_filename.c_str());

    if (!fin)
    {
        std::cerr << "Error: Could not load the centroids." << std::endl;
        exit(1);
    }

    std::vector<std::string> lines;

    std::string line;
    while (std::getline(fin, line))
        if (!line.empty())
            lines.push_back(line);

    fin.close();

    for (size_t i = 0; i < lines.size(); ++i)
    {
        std::vector<std::string> cols;

        std::string temp;
        std::istringstream stringStream;
        stringStream.str(lines[i]);

        while (stringStream >> temp)
            cols.push_back(temp);

        Point p;
        for (size_t j = 0; j < cols.size(); ++j)
            p[j] = atof(cols[j].c_str());

        centroids.push_back(p);
    }

    std::cout << "\nLoaded " << centroids.size() << " centroids.\n";

    return centroids;
}

template <typename Point>
void check_centroids(const std::vector<Point> &centroids, const size_t nb_clusters, const size_t dim)
{
    assert(nb_clusters > 0 && dim > 0);

    if (centroids.size() != nb_clusters)
    {
        std::cerr << "Error: The number of clusters "
                  << nb_clusters
                  << " is not equal to the number of loaded elements "
                  << centroids.size() << ".\n";
        exit(1);
    }

    if (centroids[0].size() != dim)
    {
        std::cerr << "Error: The number of dimensions "
                  << dim
                  << " is not equal to the dimensionality (" << centroids[0].size()
                  << ") of loaded element with index 0.\n";
        exit(1);
    }
}

// Schwefel
FIT_MAP(Schwefel){
    public :
        bool dead(){return false;
}

template <typename Indiv>
void eval(Indiv &ind)
{
    float f = 0.0f;

    for (size_t i = 0; i < ind.size(); ++i)
    {
        float internal_sum = 0.0f;

        for (size_t j = 0; j < i; ++j)
            internal_sum += ind.data(j);

        f += (internal_sum * internal_sum);
    }

    this->_value = -f;

    std::vector<float> data(Params::ea::feature_dimensionality, 0.0f);
    for (size_t i = 0; i < Params::ea::feature_dimensionality; ++i)
        data[i] = ind.data(i);

    this->set_desc(data);
}
}
;

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        std::cerr << "\nError: please provide the centroids filename.\n";
        exit(1);
    }

    Params::ea::centroids = load_centroids<point_t>(argv[1]);
    check_centroids<point_t>(Params::ea::centroids, Params::ea::number_of_clusters, Params::ea::feature_dimensionality);

    typedef Schwefel<Params> fit_t;
    typedef eval::Parallel<Params> eval_t;
    typedef gen::EvoFloat<Params::ea::dimensionality, Params> gen_t;
    typedef phen::Parameters<gen_t, fit_t, Params> phen_t;
    typedef modif::Dummy<> modifier_t;
    typedef boost::fusion::vector<stat::Map<phen_t, Params>, stat::ArchiveFit<phen_t, Params>> stat_t;

#if defined(VARIATIONGC)
    typedef ea::CVTGlobalCorrelation<phen_t, eval_t, stat_t, modifier_t, Params> ea_t;
#elif defined(VARIATIONSBX)
    typedef ea::CVTMapElites<phen_t, eval_t, stat_t, modifier_t, Params> ea_t;
#else
    typedef ea::CVTESLine<phen_t, eval_t, stat_t, modifier_t, Params> ea_t;
#endif

    ea_t ea;

    std::cout << "start run" << std::endl;
    run_ea(argc, argv, ea);
    std::cout << "end run" << std::endl;
}
