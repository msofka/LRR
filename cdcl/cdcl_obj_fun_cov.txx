#include "cdcl_obj_fun_cov.h"
#include "cdcl_utils.h"

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_transpose.h>
#include <vnl/algo/vnl_determinant.h>

#include <vcl_algorithm.h>
#include <vcl_iostream.h>

//:
// \file
// \brief  Objective function with gradients for covariance estimation.
// \author Michal Sofka
// \date   May 2006

template < unsigned int dimen, unsigned int dof >
cdcl_obj_fun_cov<dimen, dof>::cdcl_obj_fun_cov( vcl_vector< match_sptr_type > const &  matches,
                                                trans_sptr_type               const &  trans )
  : vnl_cost_function( dof*(dof+1)/2 ),
    matches_( matches ),
    trans_( trans )
{
  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< double, dimen, dimen > const &  Jp = trans_->jacobian_wrt_loc();

  // precompute mapped moving point covariances
  JCpJ_.clear();
  JCpJ_.reserve( matches_.size() );
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;

    // mapped covariance of a moving point
    vnl_matrix_fixed< double, dimen, dimen > covariance = Jp * p->covariance_ * Jp.transpose();
    
    JCpJ_.push_back( covariance );
  }
}


// Compute objective function value f and gradient g at point x.
template < unsigned int dimen, unsigned int dof >
void cdcl_obj_fun_cov<dimen, dof>::compute( const vnl_vector< double > &  x,
                                            double*                       f,
                                            vnl_vector< double >*         g )
{
  // assign parameterization into upper triangular matrix
  vnl_matrix_fixed< double, dof, dof >  upper( 0.0 );
  unsigned int ind = 0;
  for( unsigned int d1 = 0; d1 < dof; ++d1 )
    for( unsigned int d2 = d1; d2 < dof; ++d2 ) {
      upper( d1, d2 ) = x( ind );
      ++ind;
    }

  // reconstruct the covariance from the cholesky decomposition of the covariance
  trans_->set_covariance( vnl_transpose( upper )*upper );

  double num_of_all_matches = 0.0;
  vnl_matrix_fixed< double, dof, dof >  kl( 0.0 );
  if( f ) *f = 0.0;
  if( g ) g->fill( 0.0 );
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;

    // Jacobian w.r.t. transformation parameters
    vnl_matrix_fixed< double, dimen, dof > const &  Jth = trans_->jacobian_wrt_par( p->location_ );

    // transfer error and mapped covariance of a moving point
    vnl_matrix_fixed< double, dimen, dimen >  Cij_p = JCpJ_[i] + Jth * trans_->get_covariance() * Jth.transpose();
     
    num_of_all_matches += matches_[i]->to_.size();
    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches_[i]->to_[j];

      vnl_matrix_fixed< double, dimen, dimen >  Cij = Cij_p + q->covariance_;
      vnl_matrix_fixed< double, dimen, dimen >  Cij_1 = vnl_inverse( Cij );

      // matching error
      vnl_vector_fixed< double, dimen >  err_vec = q->location_ - trans_->map_loc( p->location_ );
      vnl_matrix_fixed< double, dimen, 1 >  err;
      err.set_column( 0, err_vec );

      vnl_matrix_fixed< double, 1, 1 >  residual = err.transpose() * Cij_1 * err;

      double w = matches_[i]->w_[j];

      // update objective function value
      if( f ) {
        *f += w * ( rho( residual( 0, 0 ) ) + vcl_log( vnl_determinant( Cij ) ) );
      }

      // update objective function gradient
      if( g ) {
        vnl_matrix_fixed< double, dof, 1 >  temp = Jth.transpose() * Cij_1 * err;
        kl = kl + w * ( -rho_p( residual( 0, 0 ) ) * temp * temp.transpose() + Jth.transpose() * Cij_1 * Jth );
      }
    }
  }

  if( g ) {
    // single entry matrix
    ind = 0;
    vnl_matrix_fixed< double, dof, dof >  se( 0.0 );
    for( unsigned int d1 = 0; d1 < dof; ++d1 )
      for( unsigned int d2 = d1; d2 < dof; ++d2 ) {
        se( d1, d2 ) = 1.0;
        // dCdi = d Cth/d theta_i
        vnl_matrix_fixed< double, dof, dof >  dCdth = upper.transpose() * se + se.transpose() * upper;
        (*g)[ind] = dot_product( kl.as_ref(),  dCdth.as_ref() );

        se( d1, d2 ) = 0.0;
        ++ind;
      }
  }


  if( f ) {
    assert( matches_.size() > 0 );
    *f /= num_of_all_matches;
  }

  if( g ) {
    assert( matches_.size() > 0 );
    *g *= (1.0/num_of_all_matches);
  }

  return;
}


#define CDCL_OBJ_FUN_COV_INSTANTIATE( dimen, dof ) \
template class cdcl_obj_fun_cov< dimen, dof >;
