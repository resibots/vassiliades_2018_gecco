// Vassilis Vassiliades - Inria, Nancy - 2018

#ifndef STAT_MAP_HPP_
#define STAT_MAP_HPP_

#include <boost/multi_array.hpp>
#include <limits>
#include <numeric>
#include <sferes/stat/stat.hpp>

#if !defined(VARIATIONGC) && !defined(VARIATIONSBX)
#include "archive_elem_es_line.hpp"
#endif

namespace sferes
{
namespace stat
{
SFERES_STAT(Map, Stat)
{
public:
#if !defined(VARIATIONGC) && !defined(VARIATIONSBX)
  typedef boost::shared_ptr<archive_elem_es_line<Params, Phen>> phen_t;
#else
  typedef boost::shared_ptr<Phen> phen_t;
#endif
  typedef boost::shared_ptr<Phen> indiv_t;
  typedef boost::array<float, Params::ea::feature_dimensionality> point_t;

  Map() {}

  template <typename E>
  void refresh(const E &ea)
  {
    this->_create_log_file(ea, "progress_archive.dat");

    _write_progress(ea, *this->_log_file);

    if (ea.gen() % Params::pop::dump_period == 0 || ea.gen() == Params::pop::nb_gen - 1)
    {
      _write_archive(ea.archive(), std::string("archive_"), ea);
    }
  }

protected:
  boost::shared_ptr<std::ofstream> _log_file_progress_denoised;

  template <typename E>
  void _create_custom_log_file(const E &ea, const std::string &name,
                               boost::shared_ptr<std::ofstream> &log_file)
  {
    if (!log_file && ea.dump_enabled())
    {
      std::string log = ea.res_dir() + "/" + name;
      log_file =
          boost::shared_ptr<std::ofstream>(new std::ofstream(log.c_str()));
    }
  }

  template <typename EA>
  void _write_archive(const std::vector<phen_t> &archive,
                      const std::string &prefix, const EA &ea) const
  {
    std::cout << "writing..." << prefix << ea.gen() << std::endl;
    std::string fname = ea.res_dir() + "/" + prefix +
                        boost::lexical_cast<std::string>(ea.gen()) +
                        std::string(".dat");

    std::ofstream ofs(fname.c_str());

    for (size_t i = 0; i < archive.size(); ++i)
    {
      if (archive[i])
      {
        indiv_t indiv;

#if !defined(VARIATIONGC) && !defined(VARIATIONSBX)
        indiv = archive[i]->indiv;
#else
        indiv = archive[i];
#endif
        // Write the index in archive
        ofs << i << " ";

        // Write the descriptor
        std::vector<float> descriptor = indiv->fit().desc();
        std::copy(descriptor.begin(), descriptor.end(), std::ostream_iterator<float>(ofs, " "));

        // Write the genotype
        std::vector<float> genotype = indiv->data();
        std::copy(genotype.begin(), genotype.end(), std::ostream_iterator<float>(ofs, " "));

        // Write the fitness
        ofs << " " << indiv->fit().value() << "  ";

        ofs << std::endl;
      }
    }
  }

  template <typename EA>
  void _write_progress(const EA &ea, std::ofstream &ofs) const
  {
    double archive_min = std::numeric_limits<double>::max();
    double archive_max = std::numeric_limits<double>::lowest();
    double archive_mean = 0.0;
    size_t archive_size = 0;

    std::vector<phen_t> archive = ea.archive();

    if (archive.size() == 0)
      return;

    for (size_t i = 0; i < archive.size(); ++i)
    {
      if (archive[i])
      {
        indiv_t indiv;

#if !defined(VARIATIONGC) && !defined(VARIATIONSBX)
        indiv = archive[i]->indiv;
#else
        indiv = archive[i];
#endif

        archive_size++;

        archive_mean += indiv->fit().value();

        if (indiv->fit().value() < archive_min)
          archive_min = indiv->fit().value();

        if (indiv->fit().value() > archive_max)
          archive_max = indiv->fit().value();
      }
    }

    // Divide by archive_size to calculate the mean
    archive_mean /= archive_size;

    ofs << ea.gen() << " " << ea.nb_evals() << " " << archive_size << " "
        << archive_min << " " << archive_mean << " " << archive_max
        << std::endl;
  }
};
}
}

#endif
