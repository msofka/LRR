#ifndef cdcl_trans_similarity2d_h_
#define cdcl_trans_similarity2d_h_

#include <vcl_iostream.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

#include "cdcl_trans.txx" // explicit instantiation

//:
// \file
// \brief  2D similarity transformation.
// \author Michal Sofka
// \date   May 2006

class cdcl_trans_similarity2d : public cdcl_trans< 2, 4 >
{
public:
  const static unsigned int dim = 2;
  const static unsigned int dof = 4;

  typedef vbl_smart_ptr< cdcl_trans_similarity2d >  sptr;
  typedef cdcl_trans< 2, 4 >                        superclass;
  typedef superclass::coord_type                    coord_type;

  // Default construct similarity transformation, default is identity.
  cdcl_trans_similarity2d();

  // Copy construct similarity transformation.
  cdcl_trans_similarity2d( cdcl_trans_similarity2d const &  similarity2d );

  // Construct given scale, angle, and translation.
  cdcl_trans_similarity2d( coord_type                           const &  s,
                           coord_type                           const &  ang,
                           vnl_vector_fixed< coord_type, dim >  const &  t );

  // Dynamically create similarity transformation object.
  virtual cdcl_trans< dim, dof >::sptr
  create_new() { return new cdcl_trans_similarity2d; }

  // Dynamically create copy of the current transformation object.
  virtual cdcl_trans< dim, dof >::sptr
  clone() { return new cdcl_trans_similarity2d( *this ); }

  // Return inverse of the current transform.
  virtual cdcl_trans< dim, dof >::sptr
  inverse();

  // Convert transformation to parameterization suitable for optimizer.
  virtual vnl_vector< coord_type >
  get_parameterization();

  // Convert transformation from parameterization suitable for optimizer.
  virtual void set_parameterization( vnl_vector< coord_type > const &  param );

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

  // Compute Jacobian w.r.t. inverse parameters.
  virtual vnl_matrix_fixed< coord_type, dim, dof > const &
  jacobian_wrt_inv_par( vnl_vector_fixed< coord_type, dim > const &  location );

  // Compute Jacobian w.r.t. parameters of the Jacobian w.r.t. location. 
  virtual vcl_vector< vnl_matrix_fixed< coord_type, dim, dim > > const &
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

  virtual void
  print( vcl_ostream &  ostr );

  cdcl_type_macro( cdcl_trans_similarity2d );

  // Get the rotational and scaling component (R_ matrix).
  virtual vnl_matrix_fixed< coord_type, dim, dim > &
  get_A() { A_ = s_*R_; return A_; }

  // Get the translation vector.
  virtual vnl_vector_fixed< coord_type, dim > &
  get_translation() { return t_; }

private:
  // Scale.
  coord_type                                s_;

  // Angle.
  coord_type                                ang_;

  // Rotation matrix (defined by ang_).
  vnl_matrix_fixed< coord_type, dim, dim >  R_;

  // Rotation matrix * scaling.
  vnl_matrix_fixed< coord_type, dim, dim >  A_;

  // Translation.
  vnl_vector_fixed< coord_type, dim >       t_;

};




#endif

