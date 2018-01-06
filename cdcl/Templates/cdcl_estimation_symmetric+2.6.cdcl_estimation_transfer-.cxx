#include <cdcl/cdcl_estimation_symmetric.txx>
#include <cdcl/cdcl_estimation_transfer.h>

typedef cdcl_estimation_transfer< 2, 6 > estimation_type;
CDCL_ESTIMATION_SYMMETRIC_INSTANTIATE( 2, 6, estimation_type );
