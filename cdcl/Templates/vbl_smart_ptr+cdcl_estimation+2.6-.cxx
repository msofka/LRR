#include <vbl/vbl_smart_ptr.txx>
#include <cdcl/cdcl_estimation.h>

typedef cdcl_estimation< 2, 6 >  estimation_type;
VBL_SMART_PTR_INSTANTIATE( estimation_type );
