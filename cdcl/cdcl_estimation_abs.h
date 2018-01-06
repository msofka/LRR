#ifndef cdcl_estimation_abs_h_
#define cdcl_estimation_abs_h_

#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
 
#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"

//:
// \file
// \brief  Estimation handling the whole registration: 1) matching,
//         2) parameter estimate, 3) covariance estimate.
//         Abstract interface.
// \author Michal Sofka
// \date   Nov 2006


template < unsigned int dim, unsigned int dof, class feature_typeT = cdcl_feature< dim > >
class cdcl_estimation_abs : public vbl_ref_count
{
public:
  typedef vbl_smart_ptr< cdcl_estimation_abs >  sptr;

  typedef feature_typeT                             feature_type;
  typedef typename feature_type::sptr               feature_sptr_type;
  typedef cdcl_match< dim, feature_type >           match_type;
  typedef cdcl_trans< dim, dof >                    trans_type;
  typedef typename match_type::sptr                 match_sptr_type;
  typedef typename trans_type::sptr                 trans_sptr_type;
 
  cdcl_estimation_abs() {}

  // Construct estimation object given moving and fixed feature sets
  // and an initial transformation.
  cdcl_estimation_abs( vcl_vector< feature_sptr_type > const &  moving, 
                       vcl_vector< feature_sptr_type > const &  fixed,
                       trans_sptr_type                 const &  trans );

  // Construct estimation object given moving and fixed multiresolution feature sets
  // and an initial transformation.
  // These need to have been normalized and have proper covariances.
  // This way, estimation object does not need to know anything about normalization.
  cdcl_estimation_abs( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving, 
                       vcl_vector< vcl_vector< feature_sptr_type > > const &  fixed,
                       trans_sptr_type                               const &  trans,
                       vcl_vector< double >                          const &  fixed_spacing ) {}

  typedef void   (*display_callback_type) ( void *  caller, 
                                            vcl_vector< feature_sptr_type > const &  moving_,
                                            vcl_vector< feature_sptr_type > const &  fixed_,
                                            vcl_vector< match_sptr_type >   const &  matches_,
                                            trans_sptr_type                 const &  transform,
                                            unsigned int                             iteration );

  // Initialize estimation (compute initial weights and covariance).
  virtual void initialize() = 0;

  // Run estimation.
  virtual bool run( void *  caller = 0, 
                    display_callback_type  display_points = 0 ) = 0;


};



#endif
