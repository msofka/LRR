
#include "cdcl_trans_rigid3d.h"
#include <vnl/vnl_math.h>
#include <vnl/vnl_transpose.h>
#include <vnl/algo/vnl_svd.h>

//:
// \file
// \brief  3D rigid transformation.
// \author Michal Sofka
// \date   May 2006


// Default construct rigid transformation, default is identity.
cdcl_trans_rigid3d::cdcl_trans_rigid3d()
  : cdcl_trans< dim, dof >(),
    alpha_( 0.0 ),
    beta_( 0.0 ),
    gamma_( 0.0 )
{
  t_.fill( 0.0 );
  R_.fill( 0.0 );

  R_( 0, 0 ) = 1.0;
  R_( 1, 1 ) = 1.0;
  R_( 2, 2 ) = 1.0;
}


// Construct similarity transformation given scale, angle, and translation.
cdcl_trans_rigid3d::cdcl_trans_rigid3d( coord_type                          const &  alpha,
                                        coord_type                          const &  beta,
                                        coord_type                          const &  gamma,
                                        vnl_vector_fixed< coord_type, dim > const &  t )
  : cdcl_trans< dim, dof >(),
    alpha_( alpha ),
    beta_( beta ),
    gamma_( gamma ),
    t_( t )
{
  vnl_matrix_fixed< coord_type, dim, dim >  Rx, Ry, Rz;

  coord_type sina = vcl_sin( alpha_ );
  coord_type cosa = vcl_cos( alpha_ );
  Rx(0,0) = 1.0;   Rx(0,1) = 0.0;    Rx(0,2) =  0.0;
  Rx(1,0) = 0.0;   Rx(1,1) = cosa;   Rx(1,2) = -sina;
  Rx(2,0) = 0.0;   Rx(2,1) = sina;   Rx(2,2) =  cosa;

  coord_type sinb = vcl_sin( beta_ );
  coord_type cosb = vcl_cos( beta_ );

  Ry(0,0) =  cosb;   Ry(0,1) = 0.0;   Ry(0,2) = sinb;
  Ry(1,0) =  0.0;    Ry(1,1) = 1.0;   Ry(1,2) = 0.0;
  Ry(2,0) = -sinb;   Ry(2,1) = 0.0;   Ry(2,2) = cosb;

  coord_type sing = vcl_sin( gamma_ );
  coord_type cosg = vcl_cos( gamma_ );

  Rz(0,0) = cosg;   Rz(0,1) = -sing;   Rz(0,2) = 0.0;
  Rz(1,0) = sing;   Rz(1,1) =  cosg;   Rz(1,2) = 0.0;
  Rz(2,0) = 0.0;    Rz(2,1) = 0.0;     Rz(2,2) = 1.0;
 
  R_ = Rx*Ry*Rz;
}


// Copy construct rigid transformation.
cdcl_trans_rigid3d::cdcl_trans_rigid3d( cdcl_trans_rigid3d const &  rigid3d )
  : cdcl_trans< dim, dof >( rigid3d )
{
  alpha_ = rigid3d.alpha_;
  beta_  = rigid3d.beta_;
  gamma_ = rigid3d.gamma_;

  R_     = rigid3d.R_;
  t_     = rigid3d.t_;

  covar_ = rigid3d.covar_;
}


// Construct rigid transformation given rotation matrix and translation.
cdcl_trans_rigid3d::cdcl_trans_rigid3d( vnl_matrix_fixed< coord_type, dim, dim > const &  R,
                                        vnl_vector_fixed< coord_type, dim >      const &  t,
                                        vnl_vector_fixed< coord_type, dim >      const &  center_moving )
  : R_( R ), t_( t )
{
  // determine rotation angles from the rotation matrix
  beta_    =  vcl_asin( R_(0, 2) );    // Calculate Y-axis angle
  coord_type C =  vcl_cos( beta_ );
    
  if( vcl_fabs( C ) > 1e-6 ) {             // Gimball lock?
    coord_type tr_x =  R(2,2) / C;             // No, so get X-axis angle
    coord_type tr_y = -R(1,2) / C;
    alpha_      = vcl_atan2( tr_y, tr_x );
    tr_x        =  R(0,0) / C;             // Get Z-axis angle
    tr_y        = -R(0,1) / C;
    gamma_      = vcl_atan2( tr_y, tr_x );
  }
  else {                                   // Gimball lock has occurred
    alpha_      = 0;                       // Set X-axis angle to zero
    coord_type tr_x = R(1,1);                  // And calculate Z-axis angle
    coord_type tr_y = R(1,0);
    gamma_      = vcl_atan2( tr_y, tr_x );
  }

  this->center_moving_ = center_moving;
}


// Return inverse of the current transformation.
cdcl_trans< cdcl_trans_rigid3d::dim, cdcl_trans_rigid3d::dof >::sptr
cdcl_trans_rigid3d::inverse()
{
  vnl_svd< coord_type >  svdR( R_ );
  vnl_matrix< coord_type >  invR = svdR.inverse();
  return new cdcl_trans_rigid3d( invR, -invR * t_ + this->center_moving_, vnl_vector_fixed< coord_type, dim >( 0.0 ) );
}


// Convert transformation to parameterization suitable for optimizer.
vnl_vector< cdcl_trans_rigid3d::coord_type >
cdcl_trans_rigid3d::get_parameterization()
{
  // Populate vector with parameterization.
  vnl_vector< coord_type >  param( dof, 0.0 );
  param[0] = alpha_;
  param[1] = beta_;
  param[2] = gamma_;
  param[3] = t_[0];
  param[4] = t_[1];
  param[5] = t_[2];

  return param;
}


// Convert transformation from parameterization suitable for optimizer.
void 
cdcl_trans_rigid3d::set_parameterization( vnl_vector< coord_type > const &  param )
{
  assert( param.size() == dof );

//#pragma message("cdcl_trans_rigid3d::set_parameterization:  should this be set to private variable alpha_, beta_, gamma? Then angles and R_ are ambiguous. R_ means different things for small angle and global transf. Ditch one or the other?")
alpha_ = param[0];
beta_  = param[1];
gamma_ = param[2];

  //vcl_cout << "-------------------------------SETTING-------------------------------" << vcl_endl;
  //vcl_cout << param[0] << "  " << param[1] << "  " << param[2] << "  " << param[3] << "  " << param[4] << "  " << param[5] << vcl_endl;
  //vcl_cout << "-------------------------------SETTING-------------------------------" << vcl_endl;

  coord_type dalpha  = param[0];
  coord_type dbeta   = param[1];
  coord_type dgamma  = param[2];
  
  // Rotation with small angle approximation
  R_(0,0) =  1.0;      R_(0,1) = -dgamma;   R_(0,2) =  dbeta;
  R_(1,0) =  dgamma;   R_(1,1) =  1.0;      R_(1,2) = -dalpha;
  R_(2,0) = -dbeta;    R_(2,1) =  dalpha;   R_(2,2) =  1.0;

  // Translation.
  t_[0] = param[3];
  t_[1] = param[4];
  t_[2] = param[5];
}


// Map location.
vnl_vector_fixed< cdcl_trans_rigid3d::coord_type, cdcl_trans_rigid3d::dim >
cdcl_trans_rigid3d::map_loc( vnl_vector_fixed< coord_type, dim > const &  from ) const
{
  return R_ * ( from - this->center_moving_ ) + t_;
}


// Compute Jacobian w.r.t. location.
vnl_matrix_fixed< cdcl_trans_rigid3d::coord_type, cdcl_trans_rigid3d::dim, cdcl_trans_rigid3d::dim >
cdcl_trans_rigid3d::jacobian_wrt_loc() const
{
  return R_;
}


// Compute Jacobian w.r.t. parameters.
vnl_matrix_fixed< cdcl_trans_rigid3d::coord_type, cdcl_trans_rigid3d::dim, cdcl_trans_rigid3d::dof > const &
cdcl_trans_rigid3d::jacobian_wrt_par( vnl_vector_fixed< coord_type, dim > const &  location )
{
  this->jacobian_wrt_par_thread( location, this->jacobian_wrt_par_ );

  return this->jacobian_wrt_par_;
}


// Compute Jacobian w.r.t. parameters.
void
cdcl_trans_rigid3d::jacobian_wrt_par_thread( vnl_vector_fixed< coord_type, dim > const &  location,
                                             vnl_matrix_fixed< coord_type, dim, dof >  &  jacobian_wrt_par )
{
  vnl_vector_fixed< coord_type, dim >  const  location_centered = location - this->center_moving_;

  jacobian_wrt_par( 0, 0 ) =  0.0;
  jacobian_wrt_par( 1, 0 ) = -location_centered[2];
  jacobian_wrt_par( 2, 0 ) =  location_centered[1];

  jacobian_wrt_par( 0, 1 ) =  location_centered[2];
  jacobian_wrt_par( 1, 1 ) =  0.0;
  jacobian_wrt_par( 2, 1 ) = -location_centered[0];
    
  jacobian_wrt_par( 0, 2 ) = -location_centered[1];
  jacobian_wrt_par( 1, 2 ) =  location_centered[0];
  jacobian_wrt_par( 2, 2 ) =  0.0;

  jacobian_wrt_par( 0, 3 ) = 1.0;
  jacobian_wrt_par( 1, 3 ) = 0.0;
  jacobian_wrt_par( 2, 3 ) = 0.0;

  jacobian_wrt_par( 0, 4 ) = 0.0;
  jacobian_wrt_par( 1, 4 ) = 1.0;
  jacobian_wrt_par( 2, 4 ) = 0.0;

  jacobian_wrt_par( 0, 5 ) = 0.0;
  jacobian_wrt_par( 1, 5 ) = 0.0;
  jacobian_wrt_par( 2, 5 ) = 1.0;
}


// Compute Jacobian w.r.t. parameters of the Jacobian w.r.t. location. 
vcl_vector< vnl_matrix_fixed< cdcl_trans_rigid3d::coord_type, cdcl_trans_rigid3d::dim, cdcl_trans_rigid3d::dim > > const &
cdcl_trans_rigid3d::jacobian_of_Jp()
{
  jacobian_of_Jp_.clear();

  vnl_matrix_fixed< coord_type, dim, dim >  Mkl( 0.0 );

  Mkl( 1, 2 ) = -1.0;
  Mkl( 2, 1 ) =  1.0;
  jacobian_of_Jp_.push_back( Mkl );

  Mkl( 1, 2 ) =  0.0;
  Mkl( 2, 1 ) =  0.0;

  Mkl( 0, 2 ) =  1.0;
  Mkl( 2, 0 ) = -1.0;;
  jacobian_of_Jp_.push_back( Mkl );

  Mkl( 0, 2 ) =  0.0;
  Mkl( 2, 0 ) =  0.0;;

  Mkl( 0, 1 ) = -1.0;
  Mkl( 1, 0 ) =  1.0;
  jacobian_of_Jp_.push_back( Mkl );

  return jacobian_of_Jp_;
}


// Convert transformation to normalized coordinates.
void
cdcl_trans_rigid3d::normalize( coord_type                          const &  avg_rad_moving,
                               coord_type                          const &  avg_rad_fixed,
                               vnl_vector_fixed< coord_type, dim > const &  center_moving,
                               vnl_vector_fixed< coord_type, dim > const &  center_fixed )
{
  // convert rotation matrix
  R_ = avg_rad_moving / avg_rad_fixed * R_;

  // convert translation
  t_ = t_ / avg_rad_fixed + R_ * ( center_moving - this->center_moving_ ) / avg_rad_moving - center_fixed / avg_rad_fixed;

  // THIS IS POSSIBLE ONLY WHEN R_ USES SMALL ANGLE APPROXIMATION
  // WHICH IS USUALLY NOT THE CASE WHEN NORMALIZING
  //
  // normalization matrix
  //vnl_matrix_fixed< coord_type, dof, dof >  B( 0.0 );
  //
  //...

  //this->covar_ = B * this->covar_ * B.transpose();
}


// Convert transformation back from normalized coordinates.
void
cdcl_trans_rigid3d::unnormalize( coord_type                          const &  avg_rad_moving,
                                 coord_type                          const &  avg_rad_fixed,
                                 vnl_vector_fixed< coord_type, dim > const &  center_moving,
                                 vnl_vector_fixed< coord_type, dim > const &  center_fixed )
{
  // convert rotation matrix
  R_ = avg_rad_fixed/avg_rad_moving * R_;

  // convert translation
  t_ = - R_ * center_moving + avg_rad_fixed * t_ + center_fixed;


  // normalization matrix
  //vnl_matrix_fixed< coord_type, dof, dof >  B( 0.0 );
  //
  //...

  //this->covar_ = B * this->covar_ * B;  // don't need to transpose since B symmetrical

  this->center_moving_ = center_moving;
}


// Recompose with incremental rigid transform.
void cdcl_trans_rigid3d::recompose_increment( cdcl_trans< dim, dof >::sptr const &  increment )
//void cdcl_trans_rigid3d::recompose_increment( vbl_smart_ptr< cdcl_trans< dim, dof > > &  increment )
{
  //cdcl_trans_rigid3d* ptr = dynamic_cast< cdcl_trans_rigid3d* >( increment.as_pointer() );
  cdcl_trans_rigid3d::sptr  increment_rigid3d( dynamic_cast< cdcl_trans_rigid3d* >( increment.as_pointer() ) );
  vnl_matrix_fixed< coord_type, dim, dim >  dR = increment_rigid3d->get_A();
  vnl_vector_fixed< coord_type, dim >  dt = increment_rigid3d->get_translation();

  // orthogonalize dR by setting singular values to 1
  vnl_svd< coord_type > svd( dR );
  //dR = svd.U()*vnl_transpose( svd.V() );
  dR = svd.U()*svd.V().transpose();

  // recompose transformation
  R_ = R_ * dR;
  t_ = dR * t_ + dt;
}


void
cdcl_trans_rigid3d::print( vcl_ostream &  ostr )
{
  ostr << "alpha_: " << this->alpha_ << " beta_: " << this->beta_ << " gamma_: " << this->gamma_
       << " t: " << this->t_[0] << ", " << this->t_[1] << ", " << this->t_[2] << vcl_endl;
  ostr << "R_: " << this->R_ << vcl_endl;
  ostr << "from_center: " << this->center_moving_ << vcl_endl;
  ostr << "covar: " << vcl_endl << this->get_covariance();
}

