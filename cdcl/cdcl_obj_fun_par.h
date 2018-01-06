#ifndef cdcl_obj_fun_par_h_
#define cdcl_obj_fun_par_h_

#include <vbl/vbl_ref_count.h>
#include <vnl/vnl_cost_function.h>

#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"

//:
// \file
// \brief  Objective function with gradients for parameter estimation.
// \author Michal Sofka
// \date   May 2006

// Cost function for parameter estimation (covariances fixed).
template < unsigned int dimen, unsigned int dof >
class cdcl_obj_fun_par : public vnl_cost_function, public vbl_ref_count
{
public:

  typedef typename cdcl_feature< dimen >::sptr           feature_sptr_type;
  typedef typename cdcl_match< dimen >::sptr             match_sptr_type;
  typedef typename cdcl_trans< dimen, dof >::sptr        trans_sptr_type;
  typedef typename cdcl_feature< dimen >::coord_type     coord_type;

  cdcl_obj_fun_par() : vnl_cost_function( dof ) {}

  cdcl_obj_fun_par( vcl_vector< match_sptr_type > const &  matches,
                    trans_sptr_type               const &  trans );

  // Compute objective function value f and gradient g at point x.
  void compute( const vnl_vector< double > &  x,
                double*                       f,
                vnl_vector< double >*         g );

  // Compute objective function value f and inverse gradient g at point x.
  void compute_inverse_gradient( const vnl_vector< double > &  x,
                                 double*                       f,
                                 vnl_vector< double >*         g,
                                 trans_sptr_type      const &  forward );

private:
  // Current set of matches.
  vcl_vector< match_sptr_type >  matches_;

  // Current transformation estimate.
  trans_sptr_type  trans_;

  // Transfer error JCthJ_ = Jth * Cth * Jth^T for all moving point locations.
  // We can precompute this and then save time at each iteration.
  vcl_vector< vnl_matrix_fixed< coord_type, dimen, dimen > >  JCthJ_;

};


#endif
