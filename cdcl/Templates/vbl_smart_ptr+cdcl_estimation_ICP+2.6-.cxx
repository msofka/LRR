#include <vbl/vbl_smart_ptr.txx>
#include <cdcl/cdcl_estimation_ICP.h>

typedef cdcl_estimation_ICP< 2, 6 >  estimation_type;
VBL_SMART_PTR_INSTANTIATE( estimation_type );
