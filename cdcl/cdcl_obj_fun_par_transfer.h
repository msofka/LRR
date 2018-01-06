#ifndef cdcl_obj_fun_par_transfer_h_
#define cdcl_obj_fun_par_transfer_h_

#include <vbl/vbl_ref_count.h>
#include <vnl/vnl_cost_function.h>

#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"

//:
// \file
// \brief  Objective function with gradients for parameter estimation.
//         Transfer error covariance is assumed to be estimated C_J [dimen x dimen].
// \author Michal Sofka
// \date   Oct 2006

// Cost function for parameter estimation (covariances fixed).
template < unsigned int dimen, unsigned int dof >
class cdcl_obj_fun_par_transfer : public vnl_cost_function, public vbl_ref_count
{
public:

  typedef cdcl_feature< dimen >           feature_type;
  typedef typename feature_type::sptr              feature_sptr_type;
  typedef typename cdcl_match< dimen >::sptr       match_sptr_type;
  typedef typename cdcl_trans< dimen, dof >::sptr  trans_sptr_type;

  cdcl_obj_fun_par_transfer() : vnl_cost_function( dof ) {}

  cdcl_obj_fun_par_transfer( vcl_vector< match_sptr_type > const &  matches,
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

protected:
  // Current set of matches.
  vcl_vector< match_sptr_type >  matches_;

  // Current transformation estimate.
  trans_sptr_type  trans_;

};


#endif
