#include <vbl/vbl_smart_ptr.txx>
#include <cdcl/cdcl_estimation_transfer.h>

typedef cdcl_estimation_transfer< 3, 6 >  estimation_type;
VBL_SMART_PTR_INSTANTIATE( estimation_type );
