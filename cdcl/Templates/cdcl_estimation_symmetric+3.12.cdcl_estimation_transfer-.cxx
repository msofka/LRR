#include <cdcl/cdcl_estimation_symmetric.txx>
#include <cdcl/cdcl_estimation_transfer.h>

typedef cdcl_estimation_transfer< 3, 12 > estimation_type;
CDCL_ESTIMATION_SYMMETRIC_INSTANTIATE( 3, 12, estimation_type );
