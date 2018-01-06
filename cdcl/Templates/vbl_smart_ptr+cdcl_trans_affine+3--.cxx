#include <vbl/vbl_smart_ptr.txx>
#include <cdcl/cdcl_trans_affine.h>

typedef cdcl_trans_affine< 3, 12 > trans_affine_type;
VBL_SMART_PTR_INSTANTIATE( trans_affine_type );
