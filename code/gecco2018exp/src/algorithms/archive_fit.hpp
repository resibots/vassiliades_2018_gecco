// Vassilis Vassiliades - Inria, Nancy - 2018

#include <sferes/stat/stat.hpp>
namespace sferes
{
namespace stat
{
SFERES_STAT(ArchiveFit, Stat)
{
  public:
    typedef boost::shared_ptr<Phen> indiv_t;
    static constexpr int special_value_for_empty_cell = -999999;

    ArchiveFit() {}

    template <typename E>
    void refresh(const E &ea)
    {
        if (ea.gen() % Params::pop::dump_period == 0 || ea.gen() == Params::pop::nb_gen - 1)
        {
            _write_archive_fitness(std::string("archive_fit_"), ea);
        }
    }

  protected:
    template <typename EA>
    void _write_archive_fitness(const std::string &prefix, const EA &ea) const
    {
        std::cout << "writing..." << prefix << ea.gen() << std::endl;
        std::string fname = ea.res_dir() + "/" + prefix +
                            boost::lexical_cast<std::string>(ea.gen()) +
                            std::string(".dat");

        std::ofstream ofs(fname.c_str());

        for (size_t i = 0; i < ea.archive().size(); ++i)
        {
            ofs << i << " ";

            if (ea.archive()[i]) //empty cell
            {
                indiv_t indiv;

#if !defined(VARIATIONGC) && !defined(VARIATIONSBX)
                indiv = ea.archive()[i]->indiv;
#else
                indiv = ea.archive()[i];
#endif
                ofs << indiv->fit().value();
            }
            else //non-empty cell
            {
                ofs << special_value_for_empty_cell;
            }

            ofs << std::endl;
        }
    }
};
}
}