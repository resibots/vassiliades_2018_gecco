// Vassilis Vassiliades - Inria, Nancy - 2018

#ifndef ARCHIVE_ELEM_ES_LINE_HPP_
#define ARCHIVE_ELEM_ES_LINE_HPP_

#include <cmath>

template <typename Params, typename Phen>
struct archive_elem_es_line
{
    typedef boost::shared_ptr<Phen> indiv_t;

    int parent_index; //useful to know the index of the parent when selecting
    indiv_t indiv;

    double sigma_isotropic;
    double sigma_directed;

    // Should be called only at initialization phase
    archive_elem_es_line(const indiv_t &indiv) : indiv(indiv)
    {
        parent_index = -1;
        sigma_isotropic = Params::ea::sigma_isotropic;
        sigma_directed = Params::ea::sigma_directed;
    }

    // Copy Constructor
    archive_elem_es_line(const archive_elem_es_line &rhs) : parent_index(rhs.parent_index),
                                                            sigma_isotropic(rhs.sigma_isotropic),
                                                            sigma_directed(rhs.sigma_directed)
    {
        // deep copy also the indiv so that we can modify it without modifying the original
        indiv = indiv_t(new Phen());

        for (size_t i = 0; i < indiv->gen().size(); ++i)
            indiv->gen().data(i, rhs.indiv->gen().data(i));
    }
};

#endif
