#include <vbl/vbl_smart_ptr.txx>
#include <cdcl/cdcl_estimation_symmetric.h>
#include <cdcl/cdcl_estimation_transfer.h>

typedef cdcl_estimation_transfer< 2, 4 >  estimation_type;
typedef cdcl_estimation_symmetric< 2, 4, estimation_type >  estimation_symmetric_type;
VBL_SMART_PTR_INSTANTIATE( estimation_symmetric_type );
