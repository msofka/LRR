#ifndef cdcl_trans_rigid3d_h_
#define cdcl_trans_rigid3d_h_

#include <vcl_iostream.h>
#include <vcl_cstdlib.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

#include "cdcl_trans.txx" // explicit instantiation

//:
// \file
// \brief  3D rigid transformation.
// \author Michal Sofka
// \date   May 2006

class cdcl_trans_rigid3d : public cdcl_trans< 3, 6 >
{
public:
  const static unsigned int dim = 3;
  const static unsigned int dof = 6;

  typedef vbl_smart_ptr< cdcl_trans_rigid3d >  sptr;
  typedef cdcl_trans< 3, 6 >                   superclass;
  typedef superclass::coord_type               coord_type;

  // Default construct rigid transformation, default is identity.
  cdcl_trans_rigid3d();

  // Copy construct rigid transformation.
  cdcl_trans_rigid3d( cdcl_trans_rigid3d const &  rigid3d );

  // Construct rigid transformation given scale, angles in radians, and translation.
  cdcl_trans_rigid3d( coord_type                          const &  alpha,
                      coord_type                          const &  beta,
                      coord_type                          const &  gamma,
                      vnl_vector_fixed< coord_type, dim > const &  t );

  // Construct rigid transformation given rotation matrix, translation, and a center.
  cdcl_trans_rigid3d( vnl_matrix_fixed< coord_type, dim, dim > const &  R,
                      vnl_vector_fixed< coord_type, dim >      const &  t,
                      vnl_vector_fixed< coord_type, dim >      const &  center_moving );

  // Dynamically create transformation object.
  virtual cdcl_trans< dim, dof >::sptr
  create_new() { return new cdcl_trans_rigid3d; }

  // Dynamically create copy of the current transformation object.
  virtual cdcl_trans< dim, dof >::sptr
  clone() { return new cdcl_trans_rigid3d( *this ); }

  // Return inverse of the current transform.
  virtual cdcl_trans< dim, dof >::sptr
  inverse();

  // Convert transformation to parameterization suitable for optimizer.
  virtual vnl_vector< coord_type >
  get_parameterization();

  // Convert transformation from parameterization suitable for optimizer.
  virtual void
  set_parameterization( vnl_vector< coord_type > const &  param );

  // Map location.
	virtual vnl_vector_fixed< coord_type, dim >
  map_loc( vnl_vector_fixed< coord_type, dim > const &  from ) const;

  // Compute Jacobian w.r.t. location.
  virtual vnl_matrix_fixed< coord_type, dim, dim >
  jacobian_wrt_loc() const;

  // Compute Jacobian w.r.t. parameters.
  virtual vnl_matrix_fixed< coord_type, dim, dof > const &
  jacobian_wrt_par( vnl_vector_fixed< coord_type, dim > const &  location );

  // Compute Jacobian w.r.t. parameters.
  // Thread safe version (does not modify the member variable jacobian_wrt_par_).
  virtual void
  jacobian_wrt_par_thread( vnl_vector_fixed< coord_type, dim > const &  location,
                           vnl_matrix_fixed< coord_type, dim, dof >  &  jacobian_wrt_par );

  // Compute Jacobian w.r.t. parameters of the Jacobian w.r.t. location. 
  vcl_vector< vnl_matrix_fixed< coord_type, dim, dim > > const &
  jacobian_of_Jp();

  // Convert transformation to normalized coordinates.
  virtual void 
  normalize( coord_type                          const &  avg_rad_moving,
             coord_type                          const &  avg_rad_fixed,
             vnl_vector_fixed< coord_type, dim > const &  center_moving,
             vnl_vector_fixed< coord_type, dim > const &  center_fixed );

  // Convert transformation back from normalized coordinates.
  virtual void
  unnormalize( coord_type                          const &  avg_rad_moving,
               coord_type                          const &  avg_rad_fixed,
               vnl_vector_fixed< coord_type, dim > const &  center_moving,
               vnl_vector_fixed< coord_type, dim > const &  center_fixed );

  // Recompose with incremental rigid transform.
  virtual void
  recompose_increment( cdcl_trans< dim, dof >::sptr const &  increment );
//virtual void cdcl_trans_rigid3d::recompose_increment( vbl_smart_ptr< cdcl_trans< dim, dof > > &  increment );

  // Get rotation matrix.
  virtual vnl_matrix_fixed< coord_type, dim, dim > &
  get_A() { return R_; }

  // Get translation vector.
  virtual vnl_vector_fixed< coord_type, dim > &
  get_translation() { return t_; }

  virtual void
  print( vcl_ostream &  ostr );

  cdcl_type_macro( cdcl_trans_rigid3d );

private:
  // Angles in radians w.r.t. axes.
  coord_type                                alpha_;
  coord_type                                beta_;
  coord_type                                gamma_;

  // Rotation matrix (defined by alpha_, beta_, gamma_).
  vnl_matrix_fixed< coord_type, dim, dim >  R_;

  // Translation.
  vnl_vector_fixed< coord_type, dim >       t_;

};




#endif

