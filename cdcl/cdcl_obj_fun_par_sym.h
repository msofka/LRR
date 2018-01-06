#ifndef cdcl_obj_fun_par_sym_h_
#define cdcl_obj_fun_par_sym_h_

#include <vnl/vnl_cost_function.h>

#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"
#include "cdcl_obj_fun_par.h"

//:
// \file
// \brief  Objective function with gradients for parameter estimation.
//         Forward and inverse directions (matches and covariances)
//         are considered separately. Parameter estimation function is provided
//         using forward and inverse directions consequently.
// \author Michal Sofka
// \date   May 2006

// Cost function for parameter estimation (covariances fixed).
template < unsigned int dimen, unsigned int dof, class obj_fun_par_type = cdcl_obj_fun_par< dimen, dof > >
class cdcl_obj_fun_par_sym : public vnl_cost_function
{
public:
  typedef typename cdcl_feature< dimen >::sptr     feature_sptr_type;
  typedef typename cdcl_match< dimen >::sptr       match_sptr_type;
  typedef typename cdcl_trans< dimen, dof >::sptr  trans_sptr_type;
  
  cdcl_obj_fun_par_sym() : vnl_cost_function( dof ) {}

  // Construct given matches, forward transform and inverse matches.
  cdcl_obj_fun_par_sym( vcl_vector< match_sptr_type > const &  matches,
                        trans_sptr_type               const &  trans,
                        vcl_vector< match_sptr_type > const &  matches_inverse,
                        trans_sptr_type               const &  trans_inverse );

  // Compute objective function value f and gradient g at point x.
  void compute( const vnl_vector< double > &  x,
                double*                       f,
                vnl_vector< double >*         g );

private:
  obj_fun_par_type  obj_fun_par_;
  obj_fun_par_type  obj_fun_par_inv_;

  // Current transformation estimate.
  trans_sptr_type  trans_;

};


#endif
