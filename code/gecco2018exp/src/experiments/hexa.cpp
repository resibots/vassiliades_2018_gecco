// Vassilis Vassiliades - Inria, Nancy - 2018

// #include <Eigen/Core>
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

#include <hexapod_dart/hexapod_dart_simu.hpp>

#ifdef GRAPHIC
#define NO_PARALLEL
#endif

using namespace sferes;
using namespace sferes::gen::evo_float;

struct Params
{
    struct ea
    {
        SFERES_CONST size_t dimensionality = 36;
        SFERES_CONST size_t number_of_clusters = 10000;
        SFERES_CONST size_t feature_dimensionality = 6;
        typedef boost::array<double, feature_dimensionality> point_t;
        static std::vector<point_t> centroids;

        // No self-adaptation for the directed part
        SFERES_CONST double tau_directed = 0.0;

#ifdef VARIATIONGC
        SFERES_CONST double alpha = 0.1;
#elif VARIATIONISO
        SFERES_CONST double sigma_isotropic = 0.1f;
        SFERES_CONST double sigma_directed = 0.0f;

        // No self-adaptation
        SFERES_CONST double tau_isotropic = 0.0;

        // No scaling by distance
        SFERES_CONST bool isotropic_scaled_by_distance = false;
        SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONISOSA
        SFERES_CONST double sigma_isotropic = 0.1f; // this acts as an initial value
        SFERES_CONST double sigma_directed = 0.0f;

        // Self-adaptation
        SFERES_CONST double tau_isotropic = 1.0 / sqrt(2.0 * dimensionality);

        // No scaling by distance
        SFERES_CONST bool isotropic_scaled_by_distance = false;
        SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONISODD
        SFERES_CONST double sigma_isotropic = 0.05f;
        SFERES_CONST double sigma_directed = 0.0f;

        // No self-adaptation
        SFERES_CONST double tau_isotropic = 0.0;

        // Scaling by distance for the isotropic
        SFERES_CONST bool isotropic_scaled_by_distance = true;
        SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONISOLINEDD
        SFERES_CONST double sigma_isotropic = 0.01f;
        SFERES_CONST double sigma_directed = 0.2f;

        // No self-adaptation
        SFERES_CONST double tau_isotropic = 0.0;

        // Scaling by distance for the directed
        SFERES_CONST bool isotropic_scaled_by_distance = false;
        SFERES_CONST bool directed_scaled_by_distance = true;
#elif VARIATIONLINE
    SFERES_CONST double sigma_isotropic = 0.0f;
    SFERES_CONST double sigma_directed = 0.2f;

    // No self-adaptation
    SFERES_CONST double tau_isotropic = 0.0;

    // No scaling by distance
    SFERES_CONST bool isotropic_scaled_by_distance = false;
    SFERES_CONST bool directed_scaled_by_distance = false;
#elif VARIATIONLINEDD
    SFERES_CONST double sigma_isotropic = 0.0f;
    SFERES_CONST double sigma_directed = 0.2f;

    // No self-adaptation
    SFERES_CONST double tau_isotropic = 0.0;

    // Scaling by distance for the directed
    SFERES_CONST bool isotropic_scaled_by_distance = false;
    SFERES_CONST bool directed_scaled_by_distance = true;
#endif
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

    struct pop
    {
        // number of initial random points
        SFERES_CONST size_t init_size = 100;
        SFERES_CONST unsigned size = 100;
        SFERES_CONST unsigned nb_gen = 5001;
        SFERES_CONST int dump_period = 2500;
    };
    struct parameters
    {
        SFERES_CONST float min = 0.0f;
        SFERES_CONST float max = 1.0f;
    };
};

typedef Params::ea::point_t point_t;

namespace global
{
std::shared_ptr<hexapod_dart::Hexapod> global_robot;
std::vector<int> brokenLegs;
};

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

void init_simu(std::string robot_file, std::vector<int> broken_legs = std::vector<int>())
{
    std::vector<hexapod_dart::HexapodDamage> damages(broken_legs.size());
    for (size_t i = 0; i < broken_legs.size(); ++i)
        damages.push_back(hexapod_dart::HexapodDamage("leg_removal", std::to_string(broken_legs[i])));
    global::global_robot = std::make_shared<hexapod_dart::Hexapod>(robot_file, damages);
}

FIT_MAP(FitAdapt)
{
  public:
    template <typename Indiv>
    void eval(Indiv & indiv)
    {
        this->_dead = false;
        _eval(indiv);
    }

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        dbg::trace trace("fit", DBG_HERE);

        ar &boost::serialization::make_nvp("_value", this->_value);
        ar &boost::serialization::make_nvp("_objs", this->_objs);
    }

    bool dead() { return _dead; }
    std::vector<double> ctrl() { return _ctrl; }

  protected:
    bool _dead;
    std::vector<double> _ctrl;

    template <typename Indiv>
    void _eval(Indiv & indiv)
    {
        // copy of controler's parameters
        _ctrl.clear();
        for (size_t i = 0; i < 36; i++)
            _ctrl.push_back(indiv.data(i));

        // launching the simulation
        auto robot = global::global_robot->clone();

        using desc_t = boost::fusion::vector<hexapod_dart::descriptors::DutyCycle>;
        using safe_t =
            boost::fusion::vector<hexapod_dart::safety_measures::BodyColliding,
                                  hexapod_dart::safety_measures::MaxHeight,
                                  hexapod_dart::safety_measures::TurnOver>;

        hexapod_dart::HexapodDARTSimu<hexapod_dart::safety<safe_t>,
                                      hexapod_dart::desc<desc_t>>
            simu(_ctrl, robot);

        simu.run(5);

        this->_value = simu.covered_distance();

        std::vector<float> desc(Params::ea::feature_dimensionality, 0.0f);

        if (this->_value < -1000)
        {
            this->_dead = true;
            this->_value = -1000;
        }
        else
        {

            std::vector<double> v;
            simu.get_descriptor<hexapod_dart::descriptors::DutyCycle>(v);

            desc[0] = v[0];
            desc[1] = v[1];
            desc[2] = v[2];
            desc[3] = v[3];
            desc[4] = v[4];
            desc[5] = v[5];
        }

        this->set_desc(desc);
    }
};

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        std::cerr << "\nError: please provide the centroids filename.\n";
        exit(1);
    }

    Params::ea::centroids = load_centroids<point_t>(argv[1]);
    check_centroids<point_t>(Params::ea::centroids, Params::ea::number_of_clusters, Params::ea::feature_dimensionality);

    typedef FitAdapt<Params> fit_t;
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

    std::cout << "init SIMU" << std::endl;

    const char *env_p = std::getenv("RESIBOTS_DIR");

    if (env_p) // if the environment variable exists
        init_simu(std::string(env_p) + "/share/hexapod_models/URDF/pexod.urdf",
                  global::brokenLegs);
    else // if it does not exist, we might be running this on the cluster
        init_simu(
            "/nfs/hal01/vvassili/ResiBots/share/hexapod_models/URDF/pexod.urdf",
            global::brokenLegs);

    std::cout << "start run" << std::endl;
    run_ea(argc, argv, ea);
    std::cout << "end run" << std::endl;

    return 0;
}
