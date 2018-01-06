#include <vbl/vbl_smart_ptr.txx>
#include <cdcl/cdcl_trans.h>

typedef cdcl_trans<3,12>  trans_affine3d_type;
VBL_SMART_PTR_INSTANTIATE( trans_affine3d_type );
