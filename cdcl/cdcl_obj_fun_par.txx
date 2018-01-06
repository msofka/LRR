#include "cdcl_obj_fun_par.h"
#include "cdcl_utils.h"

#include <vnl/vnl_trace.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_transpose.h>

#include <vnl/algo/vnl_determinant.h>

#include <vcl_vector.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>

//:
// \file
// \brief  Objective function with gradients for parameter estimation.
// \author Michal Sofka
// \date   May 2006

template < unsigned int dimen, unsigned int dof >
cdcl_obj_fun_par<dimen, dof>::cdcl_obj_fun_par( vcl_vector< match_sptr_type > const &  matches,
                                                trans_sptr_type               const &  trans )
  : vnl_cost_function( dof ),
    matches_( matches ),
    trans_( trans )
{
  // precompute transfer error covariance
  JCthJ_.clear();
  JCthJ_.reserve( matches_.size() );
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;

    // Jacobian w.r.t. transformation parameters
    vnl_matrix_fixed< coord_type, dimen, dof > const &  Jth = trans_->jacobian_wrt_par( p->location_ );

    // transfer error a moving point
    vnl_matrix_fixed< coord_type, dimen, dimen >  C = Jth * trans_->get_covariance() * Jth.transpose();
    JCthJ_.push_back( C );
  }
}


// Compute objective function value f and gradient g at point x.
template < unsigned int dimen, unsigned int dof >
void cdcl_obj_fun_par<dimen, dof>::compute( const vnl_vector< double > &  x,
                                            double*                       f,
                                            vnl_vector< double >*         g )
{
  // update transformation
  trans_->set_parameterization( x );

  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< coord_type, dimen, dimen > const &  Jp = trans_->jacobian_wrt_loc();

  // the following used for correcting terms coming from derivative of Jp w.r.t. parameters
  //
  // only dof-dimen elements need to be corrected because the derivative of Jp w.r.t. translation is zero
  // therefore Mkls should have size dof-dimen
  vcl_vector< vnl_matrix_fixed< coord_type, dimen, dimen > > const &  Mkls = trans_->jacobian_of_Jp();
  assert( Mkls.size() == dof-dimen );


  if( f ) *f = 0.0;
  if( g ) g->fill( 0.0 );
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;

    // Jacobian w.r.t. transformation parameters
    vnl_matrix_fixed< coord_type, dimen, dof > const &  Jth = trans_->jacobian_wrt_par( p->location_ );

    // transfer error and mapped covariance of a moving point
    vnl_matrix_fixed< coord_type, dimen, dimen >  Cij_p = JCthJ_[i] + 
                                                      Jp * p->covariance_ * Jp.transpose();
      
    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches_[i]->to_[j];

      // covariance of a correspondence match
      vnl_matrix_fixed< coord_type, dimen, dimen >  Cij = Cij_p + q->covariance_;
      vnl_matrix_fixed< coord_type, dimen, dimen >  Cij_1 = vnl_inverse( Cij );

      // matching error
      vnl_vector_fixed< coord_type, dimen >  const &  err_vec = trans_->map_loc( p->location_ ) - q->location_;
      vnl_matrix_fixed< coord_type, dimen, 1 >  err;
      err.set_column( 0, err_vec );

      vnl_matrix_fixed< coord_type, 1, 1 >  residual = err.transpose() * Cij_1 * err;
      double w = matches_[i]->w_[j];

      double rho_p_residual = rho_p( residual( 0, 0 ) );
      double rho_residual   = rho( residual( 0, 0 ) );

      // objective function or the derivative is computed only for inliers (zero contribution for outliers, because weight is zero)
      // we need to check for this because a point can become an outlier during parameter estimation
      // the correcting terms below would become invalid
      bool inlier = rho_p_residual > 0;

      // update objective function value
      if( f && inlier ) {
        *f += w * ( rho_residual + vcl_log( vnl_determinant( Cij ) ) );
      }

      // update objective function gradient
      if( g && inlier ) {
        vnl_matrix_fixed< coord_type, dimen, 1 >  eCij_1 = Cij_1 * err;
        vnl_matrix_fixed< coord_type, dimen, dimen >  dCpdth = 2.0 * p->covariance_ * Jp.transpose() * eCij_1 * eCij_1.transpose();

        vnl_matrix_fixed< coord_type, dof, 1 >  correct( 0.0 );
        vnl_matrix_fixed< coord_type, dof, 1 >  correctS( 0.0 );
        for( unsigned int f = 0; f < dof-dimen; ++f ) {
          vnl_matrix_fixed< coord_type, dimen, dimen >  Mkl = Mkls[f];

          vnl_matrix_fixed< coord_type, dimen, dimen >  dCijdth = Jp * p->covariance_ * Mkl.transpose() + Mkl * p->covariance_ * Jp.transpose();
          vnl_matrix_fixed< coord_type, 1, 1 >  dCpdth = eCij_1.transpose() * dCijdth * eCij_1;
          // correction coming form d Cij / d th  in  e^T Cij e
          correct( f, 0 ) = - dCpdth( 0, 0 );

          // correction coming from ln det Cij
          correctS( f, 0 ) = vnl_trace( Cij_1 * dCijdth );
        }

        vnl_matrix_fixed< coord_type, dof, 1 >  grad_inc = w * ( rho_p_residual * ( Jth.transpose() * Cij_1 * err * 2.0 + correct ) + correctS );
        *g += grad_inc.get_column( 0 );
      }
    }
  }

  if( f ) {
    assert( matches_.size() > 0 );
    *f /= matches_.size();
  }

  if( g ) {
    assert( matches_.size() > 0 );
    *g *= (1.0/matches_.size());
  }

  return;
}


// Compute objective function value f and inverse gradient g at point x.
template < unsigned int dimen, unsigned int dof >
void cdcl_obj_fun_par<dimen, dof>::compute_inverse_gradient( const vnl_vector< double > &  x,
                                                             double*                       f,
                                                             vnl_vector< double >*         g,
                                                             trans_sptr_type      const &  forward )
{
  // update transformation
  trans_->set_parameterization( x );

  // Jacobian w.r.t. p (moving point)
  vnl_matrix_fixed< coord_type, dimen, dimen > const &  Jp = trans_->jacobian_wrt_loc();

  // the following used for correcting terms coming from derivative of Jp w.r.t. parameters
  //
  // only dof-dimen elements need to be corrected because the derivative of Jp w.r.t. translation is zero
  // therefore Mkls should have size dof-dimen
  vcl_vector< vnl_matrix_fixed< coord_type, dimen, dimen > > const &  Mkls = trans_->jacobian_of_Jp();
  assert( Mkls.size() == dof-dimen );

  // assuming that we have linear transformations that can be written in a form y = Ax + t
  vnl_matrix_fixed< coord_type, dimen, dimen > const &  A = trans_->get_A();

  double num_of_all_matches = 0.0;
  if( f ) *f = 0.0;
  if( g ) g->fill( 0.0 );
  // go through all moving points
  for( typename vcl_vector< match_sptr_type >::size_type  i = 0; i < matches_.size(); ++i ) {
    // moving point p
    feature_sptr_type  p = matches_[i]->from_;

    // Jacobian w.r.t. transformation parameters
    vnl_vector_fixed< coord_type, dimen >   translation = forward->get_translation();
    vnl_matrix_fixed< coord_type, dimen, dof > const &  Jth_inv = trans_->jacobian_wrt_inv_par( p->location_ - translation );
      
    // transfer error and mapped covariance of a moving point
    vnl_matrix_fixed< coord_type, dimen, dimen >  Cij_p = JCthJ_[i] + Jp * p->covariance_ * Jp.transpose();
    
    num_of_all_matches += matches_[i]->to_.size();
    // go through all matches (fixed points)
    for( typename vcl_vector< match_sptr_type >::size_type  j = 0; j < matches_[i]->to_.size(); ++j ) {
      // fixed point q
      feature_sptr_type  q = matches_[i]->to_[j];

      // covariance of a correspondence match
      vnl_matrix_fixed< coord_type, dimen, dimen >  Cij   = Cij_p + q->covariance_;
      vnl_matrix_fixed< coord_type, dimen, dimen >  Cij_1 = vnl_inverse( Cij );

      // matching error
      vnl_vector_fixed< coord_type, dimen >  const &  err_vec = trans_->map_loc( p->location_ ) - q->location_;
      vnl_matrix_fixed< coord_type, dimen, 1 >  err;
      err.set_column( 0, err_vec );

      vnl_matrix_fixed< coord_type, 1, 1 >  residual = err.transpose() * Cij_1 * err;
      double w = matches_[i]->w_[j];

      double rho_p_residual = rho_p( residual( 0, 0 ) );
      double rho_residual   = rho( residual( 0, 0 ) );

      // objective function or the derivative is computed only for inliers (zero contribution for outliers, because weight is zero)
      // we need to check for this because a point can become an outlier during parameter estimation
      // the correcting terms below would become invalid
      bool inlier = rho_p_residual > 0;

      // update objective function value
      if( f && inlier ) {
        *f += w * ( rho_residual + vcl_log( vnl_determinant( Cij ) ) );
      }

      // update objective function gradient
      if( g && inlier ) {
        vnl_matrix_fixed< coord_type, dimen, 1 >  eCij_1 = Cij_1 * err;

        vnl_matrix_fixed< coord_type, dof, 1 >  correct( 0.0 );
        vnl_matrix_fixed< coord_type, dof, 1 >  correctS( 0.0 );
        for( unsigned int f = 0; f < dof-dimen; ++f ) {
          vnl_matrix_fixed< coord_type, dimen, dimen >  Mkl = Mkls[f];
          vnl_matrix_fixed< coord_type, dimen, dimen >  temp_kl = - A * Mkl * A;

          // correction coming from ln det Cij
          vnl_matrix_fixed< coord_type, dimen, dimen >  dCijdth = Jp * p->covariance_ * temp_kl.transpose() + temp_kl * p->covariance_ * Jp.transpose();
          correctS( f, 0 ) = vnl_trace( Cij_1 * dCijdth );

          // correction coming form d Cij / d th  in  e^T Cij e
          vnl_matrix_fixed< coord_type, 1, 1 >  dCpdth = eCij_1.transpose() * dCijdth * (- eCij_1);
          correct( f, 0 ) = dCpdth( 0, 0 );          
        }

        vnl_matrix_fixed< coord_type, dof, 1 >  grad_inc = w * ( rho_p_residual * ( Jth_inv.transpose() * Cij_1 * err * 2.0 + correct ) + correctS );
        *g += grad_inc.get_column( 0 );

      }
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


#define CDCL_OBJ_FUN_PAR_INSTANTIATE( dimen, dof ) \
template class cdcl_obj_fun_par< dimen, dof >;
