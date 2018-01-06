#ifndef cdcl_lbfgs_h_
#define cdcl_lbfgs_h_

//:
// \file
// \brief  Limited memory Broyden Fletcher Goldfarb Shannon minimization.
//         vnl_lbfgs minimizer enhanced with a stopping function.
//         Minimization stops when the objective function has been sufficiently decreased.
// \author Michal Sofka
// \date   Sep 2006

#include <vnl/algo/vnl_lbfgs.h>


class cdcl_lbfgs : public vnl_lbfgs
{
 public:
  //: Default constructor.
  cdcl_lbfgs() : vnl_lbfgs() {};

  //: Constructor. f is the cost function to be minimized.
  cdcl_lbfgs( vnl_cost_function &  f ) : vnl_lbfgs( f ) {};

  //: Change in the objective function since the begining of the minimization.
  double obj_fun_change() { return vcl_abs( start_error_ - end_error_ ); }

 protected:
  //: Reporting information at each iteration.
  //  Stop minimization when the objective function has been sufficiently decreased.
  virtual bool report_iter() {
    vnl_lbfgs::report_iter();
    const double tolerance = 0.05*vcl_abs( start_error_ ); // 5 percent
    vcl_cout << "START_ERROR: " << start_error_ << " END_ERROR_ " << end_error_ << " TOL: " << vcl_abs( start_error_ - end_error_ ) << vcl_endl;
    if( this->obj_fun_change() > tolerance ) return true;
    return false;
  }
  
};

#endif
