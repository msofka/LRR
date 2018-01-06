#include "cdcl_obj_fun_par_sym.h"
#include "cdcl_utils.h"

//:
// \file
// \brief  Objective function with gradients for parameter estimation.
//         Forward and inverse directions (matches and covariances)
//         are considered separately. Parameter estimation function is provided
//         using forward and inverse directions consequently.
// \author Michal Sofka
// \date   May 2006

// Construct given matches, forward transform and inverse matches.
template < unsigned int dimen, unsigned int dof, class obj_fun_par_type >
cdcl_obj_fun_par_sym<dimen, dof, obj_fun_par_type>::cdcl_obj_fun_par_sym( vcl_vector< match_sptr_type > const &  matches,
                                                        trans_sptr_type               const &  trans,
                                                        vcl_vector< match_sptr_type > const &  matches_inverse,
                                                        trans_sptr_type               const &  trans_inverse )
  : vnl_cost_function( dof ),
    trans_( trans ),
    obj_fun_par_( matches, trans ),
    obj_fun_par_inv_( matches_inverse, trans_inverse )
{
}


// Compute objective function value f and gradient g at point x.
template < unsigned int dimen, unsigned int dof, class obj_fun_par_type >
void cdcl_obj_fun_par_sym<dimen, dof, obj_fun_par_type>::compute( const vnl_vector< double > &  x,
                                                double*                       f,
                                                vnl_vector< double >*         g )
{
  // obtain parameterization of the inverse transform
  // this is the place where the consistency of forward and inverse transform is ensured
  trans_->set_parameterization( x );
  trans_sptr_type  trans_inverse = trans_->inverse();
  vnl_vector< double >  x_inverse = trans_inverse->get_parameterization();

  double f_forward = 0.0, f_inverse = 0.0;
  vnl_vector< double >  g_forward( x.size(), 0.0 ), g_inverse( x_inverse.size(), 0.0 );

  // Evaluate the objective function.
  obj_fun_par_.compute( x, &f_forward, &g_forward );
  obj_fun_par_inv_.compute_inverse_gradient( x_inverse, &f_inverse, &g_inverse, trans_ );

  if( f ) *f = f_forward + f_inverse;
  if( g ) *g = g_forward + g_inverse;

  return;
}


#define CDCL_OBJ_FUN_PAR_SYM_INSTANTIATE( dimen, dof, obj_fun_par_type ) \
template class cdcl_obj_fun_par_sym< dimen, dof, obj_fun_par_type >;
