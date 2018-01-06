#ifndef cdcl_obj_fun_cov_transfer_h_
#define cdcl_obj_fun_cov_transfer_h_

#include <vbl/vbl_ref_count.h>
#include <vnl/vnl_cost_function.h>

#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"

//:
// \file
// \brief  Objective function with gradients for covariance estimation.
//         Transfer error covariance is assumed to be estimated C_J [dimen x dimen].
// \author Michal Sofka
// \date   May 2006

// Cost function for covariance estimation (parameters fixed).
template < unsigned int dimen, unsigned int dof >
class cdcl_obj_fun_cov_transfer : public vnl_cost_function, public vbl_ref_count
{
public:
  typedef typename cdcl_feature< dimen >::sptr     feature_sptr_type;
  typedef typename cdcl_match< dimen >::sptr       match_sptr_type;
  typedef typename cdcl_trans< dimen, dof >::sptr  trans_sptr_type;

  // number of independent parameters in the covariance matrix
  const static unsigned int cov_dof = dimen*(dimen+1)/2;

  cdcl_obj_fun_cov_transfer() : vnl_cost_function( cov_dof ) {}

  cdcl_obj_fun_cov_transfer( vcl_vector< match_sptr_type > const &  matches,
                             trans_sptr_type               const &  trans );

  // Compute objective function value f and gradient g at point x.
  void compute( const vnl_vector< double > &  x,
                double*                       f,
                vnl_vector< double >*         g );

protected:
  // Current set of matches.
  vcl_vector< match_sptr_type >  matches_;

  // Current transformation estimate.
  trans_sptr_type  trans_;

  // Mapped moving point covariance JCpJ_ = Jp * Cp * Jp^T for all moving point covariance.
  // We can precompute this and then save time at each iteration.
  vcl_vector< vnl_matrix_fixed< double, dimen, dimen > >  JCpJ_;

};


#endif
