//| Copyright Inria May 2017
//| This project has received funding from the European Research Council (ERC) under
//| the European Union's Horizon 2020 research and innovation programme (grant
//| agreement No 637972) - see http://www.resibots.eu
//|
//| This software is a computer library whose purpose is to optimize continuous,
//| black-box functions. It mainly implements Gaussian processes and Bayesian
//| optimization.
//| Main repository: http://github.com/resibots/limbo
//| Documentation: http://www.resibots.eu/limbo
//|
//| This software is governed by the CeCILL-C license under French law and
//| abiding by the rules of distribution of free software.  You can  use,
//| modify and/ or redistribute the software under the terms of the CeCILL-C
//| license as circulated by CEA, CNRS and INRIA at the following URL
//| "http://www.cecill.info".
//|
//| As a counterpart to the access to the source code and  rights to copy,
//| modify and redistribute granted by the license, users are provided only
//| with a limited warranty  and the software's author,  the holder of the
//| economic rights,  and the successive licensors  have only  limited
//| liability.
//|
//| In this respect, the user's attention is drawn to the risks associated
//| with loading,  using,  modifying and/or developing or reproducing the
//| software by the user in light of its specific status of free software,
//| that may mean  that it is complicated to manipulate,  and  that  also
//| therefore means  that it is reserved for developers  and  experienced
//| professionals having in-depth computer knowledge. Users are therefore
//| encouraged to load and test the software's suitability as regards their
//| requirements in conditions enabling the security of their systems and/or
//| data to be ensured and,  more generally, to use and operate it in the
//| same conditions as regards security.
//|
//| The fact that you are presently reading this means that you have had
//| knowledge of the CeCILL-C license and that you accept its terms.
//|

#ifndef RANDOM_GENERATOR_HPP
#define RANDOM_GENERATOR_HPP

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <list>
#include <stdlib.h>
#include <random>
#include <utility>
#include <mutex>
#include <src/external/rand_utils.hpp>

namespace sferes
{
namespace tools
{
/// @ingroup tools
/// a mt19937-based random generator (mutex-protected)
///
/// usage :
/// - RandomGenerator<dist<double>>(0.0, 1.0);
/// - double r = rgen.rand();
template <typename D>
class RandomGenerator
{
  public:
    using result_type = typename D::result_type;
    RandomGenerator(result_type a, result_type b) : _dist(a, b), _rgen(randutils::auto_seed_128{}.base()) {}
    RandomGenerator(result_type a) : _dist(a), _rgen(randutils::auto_seed_128{}.base()) {}
    result_type rand()
    {
        return _dist(_rgen);
    }

  private:
    D _dist;
    std::mt19937 _rgen;
};

/// @ingroup tools
using rdist_double_t = std::uniform_real_distribution<double>;
/// @ingroup tools
using rdist_int_t = std::uniform_int_distribution<int>;
/// @ingroup tools
using rdist_gauss_t = std::normal_distribution<>;
/// @ingroup tools
using rdist_cauchy_t = std::cauchy_distribution<>;

using rdist_exponential_t = std::exponential_distribution<>;

/// @ingroup tools
/// Double random number generator
using rgen_double_t = RandomGenerator<rdist_double_t>;

/// @ingroup tools
/// Double random number generator (gaussian)
using rgen_gauss_t = RandomGenerator<rdist_gauss_t>;

/// @ingroup tools
/// Double random number generator (cauchy)
using rgen_cauchy_t = RandomGenerator<rdist_cauchy_t>;

/// @ingroup tools
/// Double random number generator (exponential)
using rgen_exponential_t = RandomGenerator<rdist_exponential_t>;

///@ingroup tools
///integer random number generator
using rgen_int_t = RandomGenerator<rdist_int_t>;

int rand_int(int min_val, int max_val)
{
    static thread_local rgen_int_t rgen(min_val, max_val);
    return rgen.rand();
}

double rand_uniform()
{
    static thread_local rgen_double_t rgen(0.0, 1.0);
    return rgen.rand();
}

double rand_normal()
{
    static thread_local rgen_gauss_t rgen(0.0, 1.0);
    return rgen.rand();
}

double rand_cauchy()
{
    static thread_local rgen_cauchy_t rgen(0.0, 1.0);
    return rgen.rand();
}

double rand_exponential()
{
    static thread_local rgen_exponential_t rgen(1.0);
    // static thread_local rdist_exponential_t d(lambda);
    return rgen.rand();
}

/// @ingroup tools
/// random vector in [0, 1.0]
///
/// - this function is thread safe because we use a random generator for each thread
/// - we use a C++11 random number generator
// Eigen::VectorXd random_vector_uniform(int size)
// {
//     static thread_local rgen_double_t rgen(0.0, 1.0);
//     Eigen::VectorXd res(size);
//     for (int i = 0; i < size; ++i)
//         res[i] = rgen.rand();
//     return res;
// }

/// @ingroup tools
/// random vector in [0, 1.0]
///
/// - this function is thread safe because we use a random generator for each thread
/// - we use a C++11 random number generator
std::vector<double> random_vector_gaussian(const int size)
{
    static thread_local rgen_gauss_t rgen(0.0, 1.0);
    std::vector<double> res(size);
    for (int i = 0; i < size; ++i)
        res[i] = rgen.rand();
    return res;
}

/// @ingroup tools
/// random vector in [0, 1.0]
///
/// - this function is thread safe because we use a random generator for each thread
/// - we use a C++11 random number generator
std::vector<double> random_vector_cauchy(const int size)
{
    static thread_local rgen_cauchy_t rgen(0.0, 1.0);
    std::vector<double> res(size);
    for (int i = 0; i < size; ++i)
        res[i] = rgen.rand();
    return res;
}
}
}

#endif
