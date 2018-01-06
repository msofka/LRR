#ifndef cdcl_estimation_symmetric_h_
#define cdcl_estimation_symmetric_h_

#include "cdcl_estimation.h"

//:
// \file
// \brief  Symmetric estimation.
//         The idea is to have forward and backward estimation object.
//         We can then run matching and covariance estimation for forward and backward case
//         and provide parameter estimation function.
// \author Michal Sofka
// \date   May 2006

template < unsigned int dim, unsigned int dof, class estimation_type = cdcl_estimation< dim, dof > >
class cdcl_estimation_symmetric : public cdcl_estimation_abs< dim, dof >
{
public:
  typedef vbl_smart_ptr< cdcl_estimation_symmetric >  sptr;

  typedef typename estimation_type::feature_type           feature_type;
  typedef typename estimation_type::feature_sptr_type      feature_sptr_type;
  typedef typename estimation_type::match_sptr_type        match_sptr_type;
  typedef typename estimation_type::trans_sptr_type        trans_sptr_type;
  typedef typename estimation_type::display_callback_type  display_callback_type;
 
  cdcl_estimation_symmetric() {}

  // Construct estimation object given moving and fixed feature sets,
  // and an initial transformation. All movind and fixed points are used.
  cdcl_estimation_symmetric( vcl_vector< feature_sptr_type >  const &  moving, 
                             vcl_vector< feature_sptr_type >  const &  fixed,
                             trans_sptr_type                  const &  trans,
                             trans_sptr_type                  const &  trans_inverse );

  // Construct estimation object given moving and fixed feature sets,
  // an initial transformation, and moving region corners.
  cdcl_estimation_symmetric( vcl_vector< feature_sptr_type >  const &  moving, 
                             vcl_vector< feature_sptr_type >  const &  fixed,
                             trans_sptr_type                  const &  trans,
                             trans_sptr_type                  const &  trans_inverse,
                             vnl_vector_fixed< double, dim >  const &  moving_x0,
                             vnl_vector_fixed< double, dim >  const &  moving_x1,
                             vnl_vector_fixed< double, dim >  const &  fixed_x0,
                             vnl_vector_fixed< double, dim >  const &  fixed_x1 );

  // Run estimation.
  virtual bool run( void *  caller = 0, 
                    display_callback_type  display_points = 0 );

  // Estimate parameters. Return true when converged.
	virtual bool estimate_parameters();

  // Estimate covariance. Return true when converged.
	virtual bool estimate_covariance();

  // Initialize estimation (compute initial weights and covariance).
  virtual void initialize();

  // Run one EM iteration of estimation (parameters, weights, covariance, weights).
  // Return true when converged.
  virtual bool one_iteration( void *  caller, 
                              display_callback_type  display_points,
                              unsigned int  iteration );

  // Set initial transform.
  void set_transforms( trans_sptr_type const &  trans,
                       trans_sptr_type const &  trans_inverse ) { forward_->set_transform( trans );
                                                                  backward_->set_transform( trans_inverse ); };

  // Set normalize_matches_ to true if the matched points should be normalized before estimation.
  virtual void normalize_matches() { normalize_matches_ = true; };

protected:

  // Normalize moving and fixed points which the matches were formed with
  // Use supplied known centers and radii
  void cdcl_normalize_matches_known( vcl_vector< match_sptr_type >        & matches,
                                   vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_moving, 
                                   typename feature_type::coord_type                                            const & avg_radius_moving, 
                                   vnl_vector_fixed< typename feature_type::coord_type, dim >                   const & center_fixed, 
                                   typename feature_type::coord_type                                            const & avg_radius_fixed );

  // Forward estimation object, estimation of forward matches and covariances.
  vbl_smart_ptr< estimation_type >  forward_;

  // Backward estimation object, estimation of forward matches and covariances.
  vbl_smart_ptr< estimation_type >  backward_;

  // True if the matched points should be normalized before estimation
  bool normalize_matches_;

};




#endif
