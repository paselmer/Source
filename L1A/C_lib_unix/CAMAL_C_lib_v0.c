#include <stdio.h>
#include <stdlib.h>

/*

This is a library of C functions intended to speed up lidar processing code.
Written by Patrick Selmer.


NOTE: TO USE THESE FUNCTIONS IN UNIX, REMOVE "__declspec" AND CHANGE
"intptr_t" TO "ssize_t."

To compile this library...
`

*/

void rebin_into_fixed_frame(double *ff, double *af, double *or , double *new,
	ssize_t *ffstrides, ssize_t *ffdims, ssize_t *afstrides, ssize_t *afdims,
	ssize_t *orstrides, ssize_t *ordims, ssize_t *newstrides, ssize_t *newdims,
	double *mult)
{
	/* Function to rebin data (counts) into a fixed frame */

	/*
	NOTES:
	-> [1/23/18]
	Core algorithm is very simple. When af bin value is between 2 ff bin values,
	assign orig value to higher of 2 ff bins. If > 1 af bins is between bin
	values, then assign mean of orig values to higher of 2 ff bins.
	It might be desireable to make this algorithm more sophistcated in the
	future.

	-> [3/6/18]
	OBSOLETE.
	*/

	/*
	ASSUMPTIONS:
	-> The input arrays are "C_CONTIGUOUS" (row major) and aligned with the
	data type such that no padded is packed between array elements.
	-> A lot of trust is being put into the calling routine (in Python)
	to get the input data types/dimensions correct. A TERRIBLE CRASH
	COULD OCCUR IF PROPER CARE ISN'T TAKEN BY CALLING ROUINE!

	INPUTS:
	ff -> 1D array of fixed frame  [fixed frame bin altitudes (ff_bins)]
	af -> 1D array of actual frame [actual frame bins (af_bins)]
	or -> the data values corresponding to the actual frame
	2D array, [number of channels x af_bins]
	new -> the data values corresponding to the fixed frame
	2D array, [number of channels x ff_bins]
	... -> strides and dimensions on each of the above 4 arrays (use
	Python ctypes convention)
	mult -> 1D array which hold # of hits in ff bin. MUST be initialized
	to all zeros. [ff_bins]

	OUTPUTS:
	The values in the memory location of "new" are replaced. Memory
	for all input variables is allocated by outside program
	(Python).
	*/


	// Declarations...

	unsigned long i, j, k, nc = ordims[0], nb_ff = ffdims[0], nb_af = afdims[0];
	unsigned long long sz_dbl;

	// Check input data. Return if assumptions not met.

	sz_dbl = sizeof(double);
	if (sz_dbl != ffstrides[0]) {
		printf("\nThe array appears to not be contiguous.\n Returning...\n");
		return;
	}

	//printf("nb_ff = %d\n", nb_ff);
	//printf("nb_af = %d\n", nb_af);

	// Perform rebinning...

	for (k = 0; k < ordims[0]; k++) {            // <- Channels

		for (j = 0; j < ffdims[0] - 1; j++) {       // <- ff_bins

			for (i = 0; i < afdims[0]; i++) {    // <- af_bins

				if (af[i] < ff[j] && af[i] >= ff[j + 1]) {
					//printf("ff:af:ff+1,i - %f:%f:%f,%d\n", ff[j],af[i],ff[j+1],i);
					new[k*nb_ff + j] = or[k*nb_af + i] + new[k*nb_ff + j];
					//if (k == 2) { printf("%d,%d %f\n", k, i, or[k*nb_af + i]); }
					mult[j] = mult[j] + 1.0;
				}

			}                                   // -> end af_bins

			if (mult[j] > 1.0) {
				new[k*nb_ff + j] = new[k*nb_ff + j] / mult[j];
			}

			mult[j] = 0.0;

		}                                      // -> end ff_bins

	}                                          // -> end Channels

}


void rebin_into_fixed_frame_v2(double *ff, double *af, double * or , double *new,
	long long *ffstrides, long long *ffdims, long long *afstrides, long long *afdims,
	long long *orstrides, long long *ordims, long long *newstrides, long long *newdims)
{
	/* Function to rebin data (counts) into a fixed frame */

	/*
	NOTES:
	-> [3/1/18]
	Reworking the original version of this algorithm/code.

	-> [3/6/18]
	Code complete
	*/

	/*
	ASSUMPTIONS:
	-> The input arrays are "C_CONTIGUOUS" (row major) and aligned with the
	data type such that no padded is packed between array elements.
	-> A lot of trust is being put into the calling routine (in Python)
	to get the input data types/dimensions correct. A TERRIBLE CRASH
	COULD OCCUR IF PROPER CARE ISN'T TAKEN BY CALLING ROUINE!

	INPUTS:
	ff -> 1D array of fixed frame  [fixed frame bin altitudes (ff_bins)]
	af -> 1D array of actual frame [actual frame bins (af_bins)]
	or -> the data values corresponding to the actual frame
	2D array, [number of channels x af_bins]
	new -> the data values corresponding to the fixed frame
	2D array, [number of channels x ff_bins]
	... -> strides and dimensions on each of the above 4 arrays (use
	Python ctypes convention)
	mult -> 1D array which hold # of hits in ff bin. MUST be initialized
	to all zeros. [ff_bins]

	OUTPUTS:
	The values in the memory location of "new" are replaced. Memory
	for all input variables is allocated by outside program
	(Python).
	*/


	// Declarations...

	unsigned long i, j, k, nc = ordims[0], nb_ff = ffdims[0], nb_af = afdims[0];
	unsigned long since_match = 0, first_match = 0, sm_lim = 2;
	unsigned long long sz_dbl;
	double mult = 0.0;

	// Check input data. Return if assumptions not met.

	sz_dbl = sizeof(double);
	if (sz_dbl != ffstrides[0] || sz_dbl != sizeof(new[0]) || sz_dbl != sizeof(or [0])) {
		printf("\nThe array appears to not be contiguous.\n Returning...\n");
		printf("\nSIZES...\nffstrids: %ld\nnewstrides: %ld\norstrides: %ld\n", ffstrides[0], newstrides[0], orstrides[0]);
		return;
	}

	//printf("nb_ff = %d\n", nb_ff);
	//printf("nb_af = %d\n", nb_af);

	// Perform rebinning...

	// Begin core of algorithm...
	for (j = 0; j < nb_ff - 1; j++) {       // <- ff_bins

		since_match = 0;
		first_match = 0;
		mult = 0.0; // # of hits in this jth bin
		i = 0;

		while (i < nb_af && since_match < sm_lim) {    // <- af_bins


			if (first_match == 1) { since_match++; }

			if (af[i] < ff[j] && af[i] >= ff[j + 1]) {

				mult = mult + 1.0;
				first_match = 1;
				since_match = 0;

				for (k = 0; k < nc; k++) {  // -> nc loop #1
					new[k*nb_ff + j] = new[k*nb_ff + j] + or [k*nb_af + i];
				}

			}

			i++;

		}                                   // -> end af_bins

		if (first_match == 1) {
			for (k = 0; k < nc; k++) {     // -> nc loop #1
				new[k*nb_ff + j] = new[k*nb_ff + j] / mult;
			}
		}

	}                                      // -> end ff_bins


}


void map_interp_times_to_orig_frame(double *NEW, double *ORIG, double *INTERP,
	double *NAV, unsigned long *NAV_MATCH, double *delta,
	ssize_t *NEWstrides, ssize_t *NEWdims, ssize_t *ORIGstrides, ssize_t *ORIGdims,
	ssize_t *INTERPstrides, ssize_t *INTERPdims, ssize_t *NAVstrides, ssize_t *NAVdims,
	ssize_t *NAV_MATCHstrides, ssize_t *NAV_MATCHdims)
{

	/* Function designed to handle a situation in which you need to overwrite
	undersampled times in data records with their interpolated and/or
	smoothed values. Also creates a map from INTERP array to ORIG & NEW arrays.
	*/

	/*
	NOTES:
	-> [2/7/18]
	-> [3/6/18]
	   The expected delta between ORIG array elements is now soft-coded as an 
	   input from the calling program (python, presumably). It had originally
	   been hard-coded in this C function as 1.0 to reflect the update rate of
	   the nav data; however, it's conceivable this update rate will change.
	*/

	/*
	ASSUMPTIONS:
	-> The input arrays are "C_CONTIGUOUS" (row major) and aligned with the
	data type such that no padded is packed between array elements.
	-> A lot of trust is being put into the calling routine (in Python)
	to get the input data types/dimensions correct. A TERRIBLE CRASH
	COULD OCCUR IF PROPER CARE ISN'T TAKEN BY CALLING ROUINE!
	-> ORIG and INTERP arrays should start at the same value.
	ORIG[0] == INTERP[0]  should be True.
	-> INTERP array increments by equal steps.

	INPUTS:
	NEW -> 1D, double-type array filled with zeros. In this array, the mapped
	interpolated values will be stored.
	ORIG -> 1D, double-type array of the original values.
	INTERP -> 1D, double-type array of the interpolated values.
	NAV_MATCH -> 1D, contains the subscripts of the NAV array that match to the
	             ORIG array.
	delta -> scalar, the expected delta between ORIG values.
	... -> strides and dimensions on each of the above 3 arrays (use
	Python ctypes convention)

	OUTPUTS:
	The values in the memory location of "NEW" are replaced. Memory
	for all input variables is allocated by outside program
	(Python).
	*/

	// Declarations...

	unsigned long i, j, k, nNEW = NEWdims[0], nORIG = ORIGdims[0], nINTERP = INTERPdims[0];
	unsigned long nNAV_MATCH = NAV_MATCHdims[0], nNAV = NAVdims[0];
	unsigned long long sz_dbl, sz_ulong;
	double del = *delta;

	// Check input data. Return if assumptions not met.

	sz_dbl = sizeof(double);
	if (sz_dbl != NEWstrides[0]) {
		printf("\nThe double array appears to not be contiguous.\n Returning...\n");
		return;
	}
	sz_ulong = sizeof(unsigned long);
	if (sz_ulong != NAV_MATCHstrides[0]) {
		printf("\nThe ulong array appears to not be contiguous.\n Returning...\n");
		return;
	}

	i = 0;
	j = 0;

	/* This while loop populates NEW with corresponding INTERP values */
	while (i < nORIG) {

		if ((INTERP[j] >= ORIG[i]) && (INTERP[j] < ORIG[i] + del)) {

			NEW[i] = INTERP[j];
			i++;
			j++;

		}
		else if (INTERP[j] < ORIG[i]) {

			// Keep scanning forward thru INTERP until INTERP[j] increases enough.
			NEW[i] = INTERP[j];
			j++;

		}
		else if (INTERP[j] > ORIG[i]) { // So, you must be greater than ORIG by at least del

		        // Repeatedly put INTERP[j] into NEW[i] until ORIG[i] increases enough.
			NEW[i] = INTERP[j];
			i++;

		}
		else {

			/* If for some reason you see the following message printed 15 times, you'll
			probably have to write some more code. */
			for (k = 0; k < 15; k++) {
				printf("\nElse statement encountered in C function map_interp_times_to_orig_frame\n");
			}
			return;

		}

		if (j >= nINTERP) { /* The interp values must not quite reach the greatest orig value.
							In this case, just keep writing the last (highest) INTERP value
							to the NEW array. */
			j--;

			while (i < nORIG) {
				NEW[i] = INTERP[j];
				i++;
			}

		}
                //printf("%f %f %d %d\n",NEW[i],INTERP[j],i,j);

	}

	/* This while loop provides an index mapping from NEW to NAV */


	i = 0; // counter for NEW in next loop
	j = 0; // counter for NAV in next loop
	while (i < nNAV_MATCH - 1) {

		if ((NAV[j] >= NEW[i]) && (NAV[j] < NEW[i + 1])) {
			NAV_MATCH[i] = j;
			i++;
			j++;
		}
		else if (NAV[j] < NEW[i]) {
			j++;
		}
		else if (NAV[j] == NEW[i]) {
			NAV_MATCH[i] = j;
			i++;
			j++;
		}
		else if (NAV[j] > NEW[i]) { // Also implied than NAV[j] >= NEW[i+1]
			NAV_MATCH[i] = j;
			i++;
		}
		else {
			/* If for some reason you see the following message printed 15 times, you'll
			probably have to write some more code. */
			for (k = 0; k < 15; k++) {
				printf("\nElse statement encountered in C function map_interp_times_to_orig_frame\n");
				printf("\nNav data is probably not meeting assumptions of this C code. Will return.\n");
			}
			return;
		}

		if (j >= nNAV) {
                        printf("In C library map_interp_times_to_orig_frame function...\n");
                        printf("There are not enough interp'd NAV data to cover ORIG data span.\n");
                        j--;
                        NAV_MATCH[i] = j;
                        // Fill in the rest of the array by repeating last match.
                        if ( i < nNAV_MATCH -1) {
                             for ( k = i; k < nNAV_MATCH - 1; k++) { NAV_MATCH[k] = j;}
                        }
			break;
		}

	}
	
	if (i == 0) { 
		printf("\n ************** WARNING **************\n");
		printf("Hello from the map_interp* function in the C library!\n");
		printf("Looks like entire NAV array is less than first NEW element.\n");
		printf("Look into this issue, please. First thing to do would be to\n");
		printf("confirm what I just said.\n");
		printf("*****************************************\n");
	}

        NAV_MATCH[nNAV_MATCH - 1] = j;	

}
