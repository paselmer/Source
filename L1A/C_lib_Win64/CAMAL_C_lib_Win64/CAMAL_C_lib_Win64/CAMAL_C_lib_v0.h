#include <stdio.h>
#include <winternl.h>

#ifdef CAMAL_C_LIB_WIN64_EXPORTS
#define CAMAL_C_LIB_WIN64_API __declspec(dllexport) 
#else
#define CAMAL_C_LIB_WIN64_API __declspec(dllimport) 
#endif

CAMAL_C_LIB_WIN64_API void rebin_into_fixed_frame(double *ff, double *af, double * or , double *NEW,
	intptr_t *ffstrides, intptr_t *ffdims, intptr_t *afstrides, intptr_t *afdims,
	intptr_t *orstrides, intptr_t *ordims, intptr_t *newstrides, intptr_t *newdims,
	double *mult);
CAMAL_C_LIB_WIN64_API void rebin_into_fixed_frame_v2(double *ff, double *af, double * or , double *NEW,
	intptr_t *ffstrides, intptr_t *ffdims, intptr_t *afstrides, intptr_t *afdims,
	intptr_t *orstrides, intptr_t *ordims, intptr_t *newstrides, intptr_t *newdims);
CAMAL_C_LIB_WIN64_API void map_interp_times_to_orig_frame(double *NEW, double *ORIG, double *INTERP,
	double *NAV, unsigned long *NAV_MATCH, double *delta,
	intptr_t *NEWstrides, intptr_t *NEWdims, intptr_t *ORIGstrides, intptr_t *ORIGdims,
	intptr_t *INTERPstrides, intptr_t *INTERPdims, intptr_t *NAVstrides, intptr_t *NAVdims,
	intptr_t *NAV_MATCHstrides, intptr_t *NAV_MATCHdims);