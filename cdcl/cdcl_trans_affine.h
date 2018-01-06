#ifndef cdcl_trans_affine_h_
#define cdcl_trans_affine_h_

#include <vcl_iostream.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>

#include "cdcl_trans.h"

//:
// \file
// \brief  Affine transformation in N dmensions.
// \author Michal Sofka
// \date   Oct 2006

// dm ... dimension
// df ... degrees of freedom
template < unsigned int dm, unsigned int df = dm*dm + dm >
class cdcl_trans_affine : public cdcl_trans< dm, df >
{
public:
  const static unsigned int dim = dm;
  const static unsigned int dof = df;

  typedef cdcl_trans< dim, df >                         superclass;
  typedef typename superclass::coord_type               coord_type;
  typedef vbl_smart_ptr< cdcl_trans_affine< dm, df > >  sptr;

  // Default construct affine transformation, default is identity.
  cdcl_trans_affine();

  // Copy construct affine transformation.
  cdcl_trans_affine( cdcl_trans_affine const &  affine );

  // Construct given affine matrix and translation.
  cdcl_trans_affine( vnl_matrix_fixed< coord_type, dm, dm > const &  A,
                     vnl_vector_fixed< coord_type, dm >     const &  t,
                     vnl_vector_fixed< coord_type, dm >     const &  center_moving );

  // Dynamically create affine transformation object.
  virtual typename cdcl_trans< dm, df >::sptr
  create_new() { return new cdcl_trans_affine<dm,df>; }

  // Dynamically create copy of the current transformation object.
  virtual typename cdcl_trans< dm, df >::sptr
  clone() { return new cdcl_trans_affine<dm,df>( *this ); }

  // Return inverse of the current transform.
  virtual vbl_smart_ptr< cdcl_trans< dm, df > >
  inverse();

  // Convert transformation to parameterization suitable for optimizer.
  virtual vnl_vector< coord_type >
  get_parameterization();

  // Convert transformation from parameterization suitable for optimizer.
  virtual void
  set_parameterization( vnl_vector< coord_type > const &  param );

  // Map location.
  virtual vnl_vector_fixed< coord_type, dm >
  map_loc( vnl_vector_fixed< coord_type, dm > const &  from ) const;

  // Compute Jacobian w.r.t. location.
  virtual vnl_matrix_fixed< coord_type, dm, dm >
  jacobian_wrt_loc() const;

  // Compute Jacobian w.r.t. parameters.
  virtual vnl_matrix_fixed< coord_type, dm, df > const &
  jacobian_wrt_par( vnl_vector_fixed< coord_type, dm > const &  location );

  // Compute Jacobian w.r.t. parameters.
  // Thread safe version (does not modify the member variable jacobian_wrt_par_).
  virtual void
  jacobian_wrt_par_thread( vnl_vector_fixed< coord_type, dm > const &  location,
                           vnl_matrix_fixed< coord_type, dm, df > &  jacobian_wrt_par );

  // Compute Jacobian w.r.t. inverse parameters.
  virtual vnl_matrix_fixed< coord_type, dm, df > const &
  jacobian_wrt_inv_par( vnl_vector_fixed< coord_type, dm > const &  location );

  // Compute Jacobian w.r.t. parameters of the Jacobian w.r.t. location. 
  virtual vcl_vector< vnl_matrix_fixed< coord_type, dm, dm > > const &
  jacobian_of_Jp();

  // Convert transformation to normalized coordinates.
  virtual void 
  normalize( coord_type                         const &  avg_rad_moving,
             coord_type                         const &  avg_rad_fixed,
             vnl_vector_fixed< coord_type, dm > const &  center_moving,
             vnl_vector_fixed< coord_type, dm > const &  center_fixed );

  // Convert transformation back from normalized coordinates.
  virtual void
  unnormalize( coord_type                         const &  avg_rad_moving,
               coord_type                         const &  avg_rad_fixed,
               vnl_vector_fixed< coord_type, dm > const &  center_moving,
               vnl_vector_fixed< coord_type, dm > const &  center_fixed );

  // Get the affine component (A matrix).
  virtual vnl_matrix_fixed< coord_type, dm, dm > &
  get_A() { return A_; }

  // Get the translation vector.
  virtual vnl_vector_fixed< coord_type, dm > &
  get_translation() { return t_; }

  virtual void
  print( vcl_ostream &  ostr );

  cdcl_type_macro( cdcl_trans_affine );

private:
  // Affine matrix.
  vnl_matrix_fixed< coord_type, dm, dm >  A_;

  // Translation.
  vnl_vector_fixed< coord_type, dm >      t_;

};


#endif

