#ifndef cdcl_estimation_h_
#define cdcl_estimation_h_

#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>

#include <vnl/vnl_transpose.h>

#include <rsdl/rsdl_kd_tree.h>

#include "cdcl_estimation_ICP.h"
#include "cdcl_obj_fun_par.h"
#include "cdcl_obj_fun_cov.h"

//:
// \file
// \brief  Estimation handling the whole registration: 1) matching,
//         2) parameter estimate, 3) covariance estimate.
//         Feature location and parameter covariances are modeled.
//         The estimation is done with expectaction maximization and covariance
//         estimated as a free parameter. Derived from ICP since it's used
//         for initialization.
// \author Michal Sofka
// \date   May 2006

template < unsigned int dim_other, unsigned int dof_other, class estimation_type >
class cdcl_estimation_symmetric;


template < unsigned int dim, unsigned int dof >
class cdcl_estimation : public cdcl_estimation_ICP< dim, dof >
{
public:
  // Symmetric estimation object is friend of estimation object
  // in order to prevent exposing private members.
  template < unsigned int dim_other, unsigned int dof_other, class estimation_type >
  friend class cdcl_estimation_symmetric;

  typedef vbl_smart_ptr< cdcl_estimation >  sptr;

  typedef cdcl_estimation_ICP< dim, dof >                  superclass_type;
  typedef typename superclass_type::feature_sptr_type      feature_sptr_type;
  typedef typename superclass_type::match_sptr_type        match_sptr_type;
  typedef typename superclass_type::trans_type             trans_type;
  typedef typename superclass_type::trans_sptr_type        trans_sptr_type;
  typedef typename superclass_type::display_callback_type  display_callback_type;
  typedef cdcl_obj_fun_cov< dim, dof >                     obj_fun_cov_type;
  typedef cdcl_obj_fun_par< dim, dof >                     obj_fun_par_type;
 
  cdcl_estimation() {}

  // Construct estimation object given moving and fixed feature sets
  // and an initial transformation.
  cdcl_estimation( vcl_vector< feature_sptr_type > const &  moving, 
                   vcl_vector< feature_sptr_type > const &  fixed,
                   trans_sptr_type                 const &  trans );

  // Construct estimation object given moving and fixed multiresolution feature sets
  // and an initial transformation.
  // These need to have been normalized and have proper covariances.
  // This way, estimation object does not need to know anything about normalization.
  cdcl_estimation( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving, 
                   vcl_vector< vcl_vector< feature_sptr_type > > const &  fixed,
                   trans_sptr_type                               const &  trans,
                   vcl_vector< double >                          const &  fixed_spacing );

  // Run estimation.
  virtual bool run( void *  caller = 0, 
                    display_callback_type  display_points = 0 );

  // Initialize estimation (compute initial weights and covariance).
  virtual void initialize();

  // Run one EM iteration of estimation (parameters, weights, covariance, weights).
  // Return true when converged.
  virtual bool one_iteration( void *  caller, 
                              display_callback_type  display_points,
                              unsigned int  iteration );

protected:
  // Compute weights, use Mahalanobis distances.
	void compute_weights( bool resolution_switch_allowed, trans_sptr_type  trans = 0 );

  // Prune matches, keep only leave_only strongest.
  void prune_matches( unsigned int leave_only );

  // Estimate parameters.
  // Return true when converged.
  virtual bool estimate_parameters();

  // Estimate covariance of parameters.
  // Return true when converged.
  virtual bool estimate_covariance();

  // Transfer error covariance.
  // The main purpose for this function is to have a single function for computing weights
  // and another for display but with different transfer error covariances.
  inline virtual vnl_matrix_fixed< typename trans_type::coord_type, dim, dim >  transfer_covar( vnl_matrix_fixed< typename trans_type::coord_type, dim, dof > const &  Jth ) {
    return Jth * this->trans_->get_covariance() * Jth.transpose();
  }
};




#endif
