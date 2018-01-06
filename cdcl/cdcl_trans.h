#ifndef cdcl_trans_h_
#define cdcl_trans_h_

#include "cdcl_macros.h"

#include <vcl_cstdlib.h>
#include <vcl_vector.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_double_2.h>

//:
// \file
// \brief  Transformation abstract class.
// \author Michal Sofka
// \date   May 2006


template < unsigned int dim, unsigned int dof >
class cdcl_trans : public vbl_ref_count
{
public:
  typedef vbl_smart_ptr< cdcl_trans<dim,dof> >      sptr;
  typedef double                                     coord_type;
  typedef vnl_matrix_fixed< coord_type, dof, dof >  covariance_type;

  // Default constructor
  cdcl_trans();

  // Copy constructor
  cdcl_trans( cdcl_trans<dim, dof> const &  trans );

  // Dynamically create a transformation object.
  virtual sptr
  create_new() = 0;

  // Dynamically create a copy of the current transformation object.
  virtual sptr
  clone() = 0;

  // Return inverse of the current transformation, if exists.
  virtual sptr
  inverse() = 0;

  // Convert transformation to parameterization suitable for optimizer.
  virtual vnl_vector< coord_type >
  get_parameterization() = 0;

  // Convert transformation from parameterization suitable for optimizer.
  virtual void
  set_parameterization( vnl_vector< coord_type > const &  param ) = 0;

  // Map location.
	virtual vnl_vector_fixed< coord_type, dim >
  map_loc( vnl_vector_fixed< coord_type, dim > const &  from ) const = 0;

  // Compute Jacobian w.r.t. location.
  virtual vnl_matrix_fixed< coord_type, dim, dim >
  jacobian_wrt_loc() const = 0;

  // Compute Jacobian w.r.t. parameters.
  virtual vnl_matrix_fixed< coord_type, dim, dof > const &
  jacobian_wrt_par( vnl_vector_fixed< coord_type, dim > const &  location ) = 0;

  // Compute Jacobian w.r.t. parameters.
  // Thread safe version (does not modify the member variable jacobian_wrt_par_).
  virtual void
  jacobian_wrt_par_thread( vnl_vector_fixed< coord_type, dim > const  &  location,
                           vnl_matrix_fixed< coord_type, dim, dof >   &  jacobian_wrt_par ) = 0;

  // Compute Jacobian w.r.t. inverse parameters.
  virtual vnl_matrix_fixed< coord_type, dim, dof > const &
  jacobian_wrt_inv_par( vnl_vector_fixed< coord_type, dim > const &  location );

  // Compute Jacobian w.r.t. parameters of the Jacobian w.r.t. location.
  virtual vcl_vector< vnl_matrix_fixed< coord_type, dim, dim > > const &
  jacobian_of_Jp() = 0;

  // Convert transformation to normalized coordinates.
  virtual void
  normalize( coord_type                          const &  avg_rad_moving,
             coord_type                          const &  avg_rad_fixed,
             vnl_vector_fixed< coord_type, dim > const &  center_moving,
             vnl_vector_fixed< coord_type, dim > const &  center_fixed ) = 0;

  // Convert transformation back from normalized coordinates.
  virtual void
  unnormalize( coord_type                          const &  avg_rad_moving,
               coord_type                          const &  avg_rad_fixed,
               vnl_vector_fixed< coord_type, dim > const &  center_moving,
               vnl_vector_fixed< coord_type, dim > const &  center_fixed ) = 0;

  // Get covariance matrix of the transformation parameter estimate.
  virtual covariance_type const &
  get_covariance() { return covar_; }

  // Get transfer error covariance estimated in the region.
  virtual vnl_matrix_fixed< coord_type, dim, dim > const &
  get_covarianceJ() { return covarJ_; }

  // Set covariance matrix of the transformation parameter estimate.
  virtual void
  set_covariance( covariance_type const &  covar ) { covar_ = covar; }

  // Set transfer error covariance estimated in the region.
  virtual void
  set_covarianceJ( vnl_matrix_fixed< coord_type, dim, dim > const &  covarJ ) { covarJ_ = covarJ; }

  // Recompose with incremental transform.
  virtual void
  recompose_increment( typename cdcl_trans< dim, dof >::sptr const &  increment ) {};

  // Get matrix part of the transformation assuming it can be written in the form Ax + t.
  virtual vnl_matrix_fixed< coord_type, dim, dim > &
  get_A() = 0;

  // Get translation vector part of the transformation assuming it can be written in the form Ax + t.
  virtual vnl_vector_fixed< coord_type, dim > &
  get_translation() = 0;

  // Print p
  virtual void
  print( vcl_ostream& ostr ) {} ;

  cdcl_type_macro( cdcl_trans );

protected:
  // Jacobian w.r.t. parameters.
  vnl_matrix_fixed< coord_type, dim, dof >  jacobian_wrt_par_;

  // Jacobian w.r.t. inverse parameters.
  vnl_matrix_fixed< coord_type, dim, dof >  jacobian_wrt_inv_par_;

  // Jacobian w.r.t. parameters of the Jacobian w.r.t. location.
  vcl_vector< vnl_matrix_fixed< coord_type, dim, dim > >  jacobian_of_Jp_;

  // Covariance matrix of the transformation parameters.
  covariance_type  covar_;

  // Transfer error covariance matrix estimated in a region Jth covar_ Jth^T.
  vnl_matrix_fixed< coord_type, dim, dim >  covarJ_;

public:
  // Store centering (center of the coordinate system in the moving space).
  vnl_vector_fixed< coord_type, dim >  center_moving_;
};


#endif

