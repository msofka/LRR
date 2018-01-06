#ifndef cdcl_trans_txx_
#define cdcl_trans_txx_

#include "cdcl_trans.h"

//:
// \file
// \brief  Transformation abstract class.
// \author Michal Sofka
// \date   May 2006


// Default constructor
template < unsigned int dim, unsigned int dof >
cdcl_trans<dim, dof>::cdcl_trans()
  : 
    jacobian_wrt_par_( 0.0 ),
    jacobian_wrt_inv_par_( 0.0 ),
    jacobian_of_Jp_( 0 ),
    covar_( 0.0 ),
    covarJ_( 0.0 ),
    center_moving_( 0.0 )
{
}


// Copy constructor
template < unsigned int dim, unsigned int dof >
cdcl_trans<dim, dof>::cdcl_trans( cdcl_trans<dim, dof> const &  trans )
  : jacobian_wrt_par_( trans.jacobian_wrt_par_ ),
    jacobian_wrt_inv_par_( trans.jacobian_wrt_inv_par_ ),
    jacobian_of_Jp_( trans.jacobian_of_Jp_ ),
    covar_( trans.covar_ ),
    covarJ_( trans.covarJ_ ),
    center_moving_( trans.center_moving_ )
{
}


// Compute Jacobian w.r.t. inverse parameters.
template < unsigned int dim, unsigned int dof >
vnl_matrix_fixed< typename cdcl_trans<dim, dof>::coord_type, dim, dof > const &
cdcl_trans<dim, dof>::jacobian_wrt_inv_par( vnl_vector_fixed< coord_type, dim > const &  location )
{
  // All children should overwrite.
  vcl_cout << "Error: Transformation does not implement jacobian_wrt_inv_par." << vcl_endl;
  vcl_abort();
  return jacobian_wrt_inv_par_;
};



#endif
