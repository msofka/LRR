#ifndef cdcl_match_h_
#define cdcl_match_h_

#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vcl_vector.h>

#include "cdcl_feature.h"

//:
// \file
// \brief  Match for one feature. This match can be one-to-many
//         and weights of these matches are also stored.
// \author Michal Sofka
// \date   May 2006

template < unsigned int dim, class feature_typeT = cdcl_feature< dim > >
class cdcl_match : public vbl_ref_count
{
public:
  //static const unsigned int dim = feature_type::dim;
  typedef vbl_smart_ptr< cdcl_match >  sptr;
  typedef feature_typeT                  feature_type;
  typedef typename feature_type::sptr    feature_sptr_type;

  cdcl_match()
    : from_( 0 ) { to_.clear(); w_.clear(); }

  cdcl_match( feature_sptr_type               const &  from, 
              vcl_vector< feature_sptr_type > const &  to,
              vcl_vector< double >            const &  w )
    : from_( from ), to_( to ), w_( w ) 
  { }

  // ``From'' (moving) point.
  feature_sptr_type                from_;

  // Vector of ``to'' (fixed) points.
  vcl_vector< feature_sptr_type >  to_;

  // Weights corresponding to each match.
  vcl_vector< double >             w_;

};



#endif
