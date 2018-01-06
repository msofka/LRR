#include <vbl/vbl_smart_ptr.txx>
#include <cdcl/cdcl_match.h>
#include <cdcl/cdcl_feature_with_shape.h>

typedef cdcl_feature_with_shape< 3 >  feature_type;
typedef cdcl_match< 3, feature_type >  match_type;
VBL_SMART_PTR_INSTANTIATE( match_type );

