#include <vcl_utility.h>
#include <vcl_iomanip.h>
#include <vcl_iostream.h>
#include <vcl_iterator.h>
#include <vcl_cstdlib.h>
#include <vcl_algorithm.h>

#include <vul/vul_timer.h>
#include <vnl/algo/vnl_cholesky.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_trace.h>
#include <vnl/vnl_transpose.h>

#include "cdcl_estimation_transfer.h"
#include "cdcl_utils.h"
#include "cdcl_lbfgs.h"


// VTK is used to dump data for analysis in external programs.
#ifdef VTK_FOUND
#include "cdcl_utils_VTK.h"
#endif

//:
// \file
// \brief  Estimation handling the whole registration: 1) matching,
//         2) parameter estimate, 3) covariance estimate.
//         Feature location covariances and transfer error covariance (one for 
//         all points are modeled.
//         The estimation is done with expectaction maximization and covariance
//         estimated as a free parameter.
//         Inheritted from cdcl_estimation which estimates full covariance. The difference
//         is in the objective functions.
// \author Michal Sofka
// \date   Nov 2006


template < unsigned int dim, unsigned int dof >
void cdcl_estimation_transfer<dim, dof>::initialize()
{

  cdcl_estimation<dim, dof>::initialize();

  bool resolution_switch_allowed = false;

  // run estimation of the full covariance Cth for 5 iterations and then initialize transfer covariance CJ
  // since we overwrote function transfer_covar(.), calls to cdcl_estimation function will use Cth
  // calls to cdcl_estimation_transfer functions will use CJ
  for( unsigned int i = 0; i < 5; ++i ) {
    cdcl_estimation<dim, dof>::compute_weights( resolution_switch_allowed );
    cdcl_estimation<dim, dof>::estimate_covariance();
  }


//  cdcl_estimation<dim, dof>::initialize();
//
//return;

  //for( unsigned int i = 0; i < 5; ++i ) {
  //  //call to cdcl_estimation<dim, dof>::one_iteration( 0, 0, 0 ); will result in calling virtual functions of cdcl_estimation_transfer
  //  bool resolution_switch_allowed = false;
  //  //vcl_cout << "COV!!!" << vcl_endl << this->trans_->get_covariance() << vcl_endl;
  //  cdcl_estimation<dim, dof>::compute_weights( resolution_switch_allowed );
  //  cdcl_estimation<dim, dof>::estimate_covariance();
  //}
  ////vcl_cout << "COV!!!" << vcl_endl << this->trans_->get_covariance() << vcl_endl;
  
  vnl_matrix_fixed< double, dim, dim >  covarianceJ;
  
  // initial transfer error covariance is the one with maximum trace of Jth Cth Jth^T, points sampled on a circle
  // double max_trace = 0.0;
  //for( double ang = 0; ang < 2*vnl_math::pi; ang+=vnl_math::pi/10.0 ) {
  //  double radius = 1.5;
  //  double loc_x = radius * vcl_cos( ang );
  //  double loc_y = radius * vcl_sin( ang );
  //  vnl_matrix_fixed< double, dim, dof >  Jth = this->trans_->jacobian_wrt_par( vnl_vector_fixed< double, dim >( loc_x, loc_y ) );
  //  vnl_matrix_fixed< double, dim, dim >  initial_covarianceJ =  Jth * this->trans_->get_covariance() * vnl_tranpose( Jth );
  //  double trace = vnl_trace( initial_covarianceJ );
  //  if( trace > max_trace ) {
  //    trace = max_trace;
  //    covarianceJ = initial_covarianceJ;
  //  }
  //}

  // initial transfer error covariance is the one with maximum trace of Jth Cth Jth^T
  // double max_trace = 0.0;
  //for( typename vcl_vector< double >::size_type  i = 0; i < moving_[res_].size(); ++i ) {
  //  // moving point p
  //  feature_sptr_type  p = moving_[res_][i];

  //  // Jacobian w.r.t. transformation parameters
  //  vnl_matrix_fixed< double, dim, dof > const &  Jth = this->trans_->jacobian_wrt_par( p->location_ );

  //  vnl_matrix_fixed< double, dim, dim >  initial_covarianceJ =  Jth * this->trans_->get_covariance() * vnl_tranpose( Jth );
  //  double trace = vnl_trace( initial_covarianceJ );
  //  if( trace > max_trace ) {
  //    trace = max_trace;
  //    covarianceJ = initial_covarianceJ;
  //  }
  //}


//  // initial transfer error covariance is average as determined by eigenvectors and eigenvalues of Jth Cth Jth^T
//  vnl_matrix_fixed< double, dim, dim >  U( 0.0 ), W( 0.0 );
//  for( typename vcl_vector< double >::size_type  i = 0; i < this->moving_[this->res_].size(); ++i ) {
//    // moving point p
//    feature_sptr_type  p = this->moving_[this->res_][i];
//
//    // Jacobian w.r.t. transformation parameters
//    vnl_matrix_fixed< double, dim, dof > const &  Jth = this->trans_->jacobian_wrt_par( p->location_ );
//
//    vnl_matrix_fixed< double, dim, dim >  initial_covarianceJ = Jth * this->trans_->get_covariance() * Jth.transpose();
//
//    vnl_svd< double > svd( initial_covarianceJ );
//    U += svd.U();
//    W += svd.W();
////    V += svd.V();
//  }
//  U *= 1.0 / this->moving_[this->res_].size();
//  W *= 1.0 / this->moving_[this->res_].size();
//  //svd.V() = V * 1.0 / moving_[res_].size();
//  covarianceJ = 1000.0*U*W*U.transpose();
  

//  // POSSIBLE PROBLEMS: when computing maximum, eigenvectors might not end up orthogonal
//  // initial transfer error covariance is maximum as determined by eigenvectors and eigenvalues of Jth Cth Jth^T
//  vnl_matrix_fixed< double, dim, dim >  U( 0.0 ), W( 0.0 );
//  for( typename vcl_vector< double >::size_type  i = 0; i < this->moving_[this->res_].size(); ++i ) {
//    // moving point p
//    feature_sptr_type  p = this->moving_[this->res_][i];
//
//    // Jacobian w.r.t. transformation parameters
//    vnl_matrix_fixed< double, dim, dof > const &  Jth = this->trans_->jacobian_wrt_par( p->location_ );
//
//    vnl_matrix_fixed< double, dim, dim >  initial_covarianceJ = Jth * this->trans_->get_covariance() * Jth.transpose();
//
//    vnl_svd< double > svd( initial_covarianceJ );
//    U += svd.U();
//    if( W( 0, 0 ) < svd.W( 0, 0 ) ) {
//      W( 0, 0 ) = svd.W( 0, 0 );
////      U.set_column( 0, svd.U().get_column( 0 ) );
//    }
//    if( W( 1, 1 ) < svd.W( 1, 1 ) ) {
//      W( 1, 1 ) = svd.W( 1, 1 ); 
////      U.set_column( 1, svd.U().get_column( 1 ) );
//    }
//  }
//  U *= 1.0 / this->moving_[this->res_].size();
//  covarianceJ = U*W*U.transpose();


  // initial transfer error covariance is average as determined by averaging coefficients of Jth Cth Jth^T
  vnl_matrix_fixed< double, dim, dim >  initial_covarianceJ( 0.0 );
  for( typename vcl_vector< double >::size_type  i = 0; i < this->moving_[this->res_].size(); ++i ) {
    // moving point p
    feature_sptr_type  p = this->moving_[this->res_][i];

    // Jacobian w.r.t. transformation parameters
    vnl_matrix_fixed< double, dim, dof > const &  Jth = this->trans_->jacobian_wrt_par( p->location_ );

    initial_covarianceJ += Jth * this->trans_->get_covariance() * Jth.transpose();
  }
  covarianceJ = 1.0 * initial_covarianceJ * (1.0 / this->moving_[this->res_].size() );


  this->trans_->set_covarianceJ( covarianceJ );
}


// Estimate parameters.
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation_transfer<dim, dof>::estimate_parameters()
{
  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( this->normalize_matches_ ) {
    this->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    this->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
  }

  // objective function
  obj_fun_par_type  f_par( this->matches_, this->trans_ );

  vnl_vector< double >  x = this->trans_->get_parameterization();  // can afford to copy the vector
  vnl_vector< double >  x_old = x;

  // estimator
  cdcl_lbfgs  minimizer( f_par );
  minimizer.set_max_function_evals( 5 );

  //minimizer.set_f_tolerance( 0.001 );
  //minimizer.set_x_tolerance( 1e-3 );
  minimizer.set_verbose( true );
  minimizer.set_check_derivatives( false );
  minimizer.set_trace( false );
  //minimizer.default_step_length = 0.1;
  bool converged = minimizer.minimize( x );
  bool parameters_converged = ( minimizer.get_start_error() - minimizer.get_end_error() ) < 1e-5;

  if( ( !converged && minimizer.get_failure_code() < 0 ) || x.size() == 0 ) {
    vcl_cerr << "Error during minimization, failure code: " << minimizer.get_failure_code() << " . Keeping the old value." << vcl_endl;
    this->trans_->set_parameterization( x_old );
    //exit( 1 );
  }
  else {
    this->trans_->set_parameterization( x );
  }

  if( this->normalize_matches_ ) {
    // convert transformation back to unnormalized coordinates
    this->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    this->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
  }

  return parameters_converged;
}


// Estimate covariance.
// Return true when converged.
template < unsigned int dim, unsigned int dof >
bool cdcl_estimation_transfer<dim, dof>::estimate_covariance()
{
  // Normalize before estimation.
  //
  // Normalization parameters.
  vnl_vector_fixed< double, dim >  center_moving( 0.0 ), center_fixed( 0.0 );
  double avg_rad_moving( 0.0 ), avg_rad_fixed( 0.0 );

  if( this->normalize_matches_ ) {
    this->cdcl_normalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
    this->trans_->normalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
  }

  // objective function
  obj_fun_cov_type  f_cov( this->matches_, this->trans_ );

  // number of independent parameters in the covariance matrix
  const unsigned int cov_dof = dim*(dim+1)/2;

  // assign upper triangular matrix, cholesky decomposition of the covariance, to x
  vnl_cholesky  chol( this->trans_->get_covarianceJ() );
  vnl_matrix< double >  upper = chol.upper_triangle();
  vnl_vector< double >  x( cov_dof, 0.0 );
  unsigned int ind = 0;
  for( unsigned int d1 = 0; d1 < dim; ++d1 )
    for( unsigned int d2 = d1; d2 < dim; ++d2 ) {
      x( ind ) = upper( d1, d2 );
      //vcl_cout << d1 << " " << d2 << " " << upper( d1, d2 ) << vcl_endl;
      ++ind;
    }

  vnl_vector< double >  x_old = x;

  // estimator
  cdcl_lbfgs  minimizer( f_cov );
  minimizer.set_max_function_evals( 5 );
  minimizer.set_verbose( true );
  minimizer.set_check_derivatives( false );
  minimizer.set_trace( false );
  //minimizer.default_step_length = 0.1;
  bool converged = minimizer.minimize( x );
  bool covariance_converged = ( minimizer.get_start_error() - minimizer.get_end_error() ) < 1e-5;

  if( ( !converged && minimizer.get_failure_code() < 0 ) || x.size() == 0 ) {
      vcl_cerr << "Error during minimization, failure code: " << minimizer.get_failure_code() << " . Keeping the old value." << vcl_endl;
      x = x_old;
      //exit( 1 );
  }

  // assign parameterization into upper triangular matrix
  upper.fill( 0.0 );
  ind = 0;
  for( unsigned int d1 = 0; d1 < dim; ++d1 )
    for( unsigned int d2 = d1; d2 < dim; ++d2 ) {
      upper( d1, d2 ) = x( ind );
      ++ind;
    }

  // reconstruct the covariance from the cholesky decomposition of the covariance
  this->trans_->set_covarianceJ( vnl_transpose( upper )*upper );

  //if( !minimized ) 
  //  vcl_abort();
  // covariance has already been reconstructed from the cholesky decomposition in the optimizer loop


  if( this->normalize_matches_ ) {
    // convert transformation back to unnormalized coordinates
    this->trans_->unnormalize( avg_rad_moving, avg_rad_fixed, center_moving, center_fixed );
    this->cdcl_unnormalize_matches( center_moving, avg_rad_moving, center_fixed, avg_rad_fixed );
  }

  return covariance_converged;
}




#define CDCL_ESTIMATION_TRANSFER_INSTANTIATE( dim, dof ) \
template                                        \
class cdcl_estimation_transfer< dim, dof >;
