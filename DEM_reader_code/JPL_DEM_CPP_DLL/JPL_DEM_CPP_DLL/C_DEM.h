///////////////////////////////////////////////////////////////////////////////
//
//	C_DEM.h  --> Edited by Patrick Selmer for Unix compatibility
//

#ifndef	_INC_C_DEM
#define	_INC_C_DEM

// ********************************************************************************************

//#include <afx.h>
#include <math.h>

/* All code below this line should match the Windows C_DEM.h version <------------------------------- */ 

// ********************************************************************************************

#define CIRAERR_UNEXPECTED_FILESIZE								4000000226

// ********************************************************************************************

#define DEM_NUM_ROWS			39840
#define DEM_NUM_COLS			86400

#define DEM_LAT_BEG				83
#define DEM_LAT_END				-83

#define DEM_LON_BEG				-180
#define DEM_LON_END				180

#define DEM_BYTES_PER_POINT		2

#define DEM_ELEV_ERROR			9999

// ********************************************************************************************

class C_DEM
{
public:
	C_DEM();
	~C_DEM();

	DWORD Init( const char *Filename, BOOL ReadIntoMemory = FALSE );

	ULONGLONG Lat2Row( double Lat );
	ULONGLONG Lon2Col( double Lon );

	double Row2Lat( ULONGLONG Row );
	double Col2Lon( ULONGLONG Col );

	short GetElevationShort( ULONGLONG Row, ULONGLONG Col, BYTE *LOCI = NULL );

	short GetElevationShort( double Lat, double Lon, BYTE *LOCI = NULL );
	short GetElevationShort( float Lat, float Lon, BYTE *LOCI = NULL );

	float GetElevationFloat( double Lat, double Lon, BYTE *LOCI = NULL );
	float GetElevationFloat( float Lat, float Lon, BYTE *LOCI = NULL );

private:

	__int64 Round64( double Number );

	FILE *m_fpIn;

	ULONGLONG COLS_PER_ROW;
	ULONGLONG BYTES_PER_ROW;

	ULONGLONG MAX_ROW_INDEX;
	ULONGLONG MAX_COL_INDEX;

	double MIN_LAT;
	double MIN_LON;

	double MAX_LAT;
	double MAX_LON;

	double ABS_MIN_LAT;
	double ABS_MIN_LON;

	double ABS_MAX_LAT;
	double ABS_MAX_LON;

	double LAT_DENOMINATOR;
	double LON_DENOMINATOR;

	double LAT_FACTOR;
	double LON_FACTOR;

	unsigned short LOCI_BITMASK;

	BOOL m_ReadIntoMemory;

	unsigned short *m_FileArray;
};

#endif	// _INC_C_DEM

/* Between '+''s added by Patrick Selmer
++++++++++++++++++++++++++++++++++++++++++++++
*/

#ifndef _DEM_OUTPUTS
struct DEM_outputs
{
	float Elev_rtn;
	unsigned char LOCI_rtn;
};

#endif //_DEM_OUTPUTS

#define DEM_FILE_READ_ERROR 9999999

/*
++++++++++++++++++++++++++++++++++++++++++++++
*/

