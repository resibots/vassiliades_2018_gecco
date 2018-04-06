// Vassilis Vassiliades - Inria, Nancy - 2018

#include <sferes/fit/fitness.hpp>

#define FIT_MAP(Name) SFERES_FITNESS(Name, sferes::fit::FitMap)

namespace sferes
{
namespace fit
{
SFERES_FITNESS(FitMap, sferes::fit::Fitness)
{
public:
  FitMap()
      : _desc(Params::ea::feature_dimensionality, 0.0f)
  {
  }

  template <typename Indiv>
  float dist(const Indiv &ind)
  {
    assert(_desc.size() == ind.fit()._desc.size() &&
           _desc.size() == Params::ea::feature_dimensionality);
    float sum = 0.0;

    for (size_t i = 0; i < _desc.size(); ++i)
    {
      float dif = desc(i) - ind.fit().desc(i);
      sum += dif * dif;
    }

    return sqrtf(sum);
  }

  const std::vector<float> &desc() const { return _desc; }

  float desc(const size_t index) const
  {
    assert(index < _desc.size());

    return _desc[index];
  }

  void set_desc(const std::vector<float> &x)
  {
    assert(x.size() == Params::ea::feature_dimensionality);

    _desc = x;
  }

  void set_value(const float &val) { this->_value = val; }

protected:
  std::vector<float> _desc;
};
}
}
