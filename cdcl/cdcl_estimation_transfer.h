#ifndef cdcl_estimation_transfer_h_
#define cdcl_estimation_transfer_h_

#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>
 
#include <rsdl/rsdl_kd_tree.h>

#include "cdcl_estimation.h"
#include "cdcl_feature.h"
#include "cdcl_match.h"
#include "cdcl_trans.h"
#include "cdcl_obj_fun_par_transfer.h"
#include "cdcl_obj_fun_cov_transfer.h"

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

template < unsigned int dim_other, unsigned int dof_other, class estimation_type >
class cdcl_estimation_symmetric;


template < unsigned int dim, unsigned int dof >
class cdcl_estimation_transfer : public cdcl_estimation< dim, dof >
{
public:
  // Symmetric estimation object is friend of estimation object
  // in order to prevent exposing private members.
  template < unsigned int dim_other, unsigned int dof_other, class estimation_type >
  friend class cdcl_estimation_symmetric;

  typedef vbl_smart_ptr< cdcl_estimation_transfer >  sptr;

  typedef cdcl_estimation< dim, dof >                      superclass_type;
  typedef typename superclass_type::feature_sptr_type      feature_sptr_type;
  typedef typename superclass_type::match_sptr_type        match_sptr_type;
  typedef typename superclass_type::trans_type             trans_type;
  typedef typename superclass_type::trans_sptr_type        trans_sptr_type;
  typedef typename superclass_type::display_callback_type  display_callback_type;
  typedef cdcl_obj_fun_cov_transfer< dim, dof >            obj_fun_cov_type;
  typedef cdcl_obj_fun_par_transfer< dim, dof >            obj_fun_par_type;
 
  cdcl_estimation_transfer() {}

  // Construct estimation object given moving and fixed feature sets
  // and an initial transformation.
  cdcl_estimation_transfer( vcl_vector< feature_sptr_type > const &  moving, 
                            vcl_vector< feature_sptr_type > const &  fixed,
                            trans_sptr_type                 const &  trans )
    : cdcl_estimation<dim,dof>( moving, fixed, trans ) {};

  // Construct estimation object given moving and fixed multiresolution feature sets
  // and an initial transformation.
  // These need to have been normalized and have proper covariances.
  // This way, estimation object does not need to know anything about normalization.
  cdcl_estimation_transfer( vcl_vector< vcl_vector< feature_sptr_type > > const &  moving, 
                            vcl_vector< vcl_vector< feature_sptr_type > > const &  fixed,
                            trans_sptr_type                               const &  trans,
                            vcl_vector< double >                          const &  fixed_spacing )
    : cdcl_estimation<dim,dof>( moving, fixed, trans, fixed_spacing ) {};

  // Initialize estimation (compute initial weights and covariance).
  virtual void initialize();

protected:
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
    return this->trans_->get_covarianceJ();
  }

};



#endif
