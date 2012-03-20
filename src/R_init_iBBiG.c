#include <R_ext/Rdynload.h>
#include "iBBiG.h"

void R_init_iBBiG(DllInfo *info)
{
    R_NativePrimitiveArgType cc_type[11] =
        {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP,
        REALSXP, REALSXP, INTSXP, REALSXP};

    R_CMethodDef cMethods[] = {
        {"clusterCovsC", (DL_FUNC) &clusterCovsC, 11, cc_type}, 
        {NULL, NULL, 0}
    };

    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
