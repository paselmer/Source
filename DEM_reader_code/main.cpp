#include <afx.h>

#include "C:\C_DEM.h"

////////////////////////////////////////////////////////////////////////////////

void main( int argc, char *argv[] )
{
	C_DEM DEM;

	DWORD Status = DEM.Init( "C:\\JPL-CloudSat.dem", FALSE );
	if( Status == ERROR_SUCCESS )
	{
		BYTE LOCI;

		float Lat = float( 40.5592 );
		float Lon = float( -105.0781 );

		short Elev = DEM.GetElevationShort( Lat, Lon, &LOCI );

		// NOTE: the raw LOCI flag extracted directly from the DEM file is different than the LOCI flag that is returned from the read routine.
		// the raw LOCI flag is a system flag designed to save space (see DEM_description.txt for more info) and the LOCI flag returned here is the finalized science version of the flag.

		CString LOCI_Description;

		if( LOCI == 1 )
			LOCI_Description = CString( "land" );
		else if( LOCI == 2 )
			LOCI_Description = CString( "ocean" );
		else if( LOCI == 3 )
			LOCI_Description = CString( "coast" );
		else if( LOCI == 4 )
			LOCI_Description = CString( "inland waters" );
		else if( LOCI == 5 )
			LOCI_Description = CString( "inland mixed" );

		printf( "Elev = %d | LOCI = %d (%s)\n", Elev, LOCI, LOCI_Description );
	}
}

