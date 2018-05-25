///////////////////////////////////////////////////////////////////////////////
//
//	C_DEM.cpp
//
//	This class is used to access a digital elevation map (DEM) in order to
//  get the elevation at a specified latitude and longitude
//
//	Notes:
//		Returned elevation of 9999 indicates an error
//

#include "C_DEM.h"

///////////////////////////////////////////////////////////////////////////////

C_DEM::C_DEM()
{
	COLS_PER_ROW = ULONGLONG( DEM_NUM_COLS );
	BYTES_PER_ROW = ULONGLONG( DEM_NUM_COLS ) * ULONGLONG( DEM_BYTES_PER_POINT );

	MAX_ROW_INDEX = ULONGLONG( DEM_NUM_ROWS - 1 );
	MAX_COL_INDEX = ULONGLONG( DEM_NUM_COLS - 1 );

	int RowsPerDeg = DEM_NUM_ROWS / (abs( DEM_LAT_BEG ) + abs( DEM_LAT_END ));
	int ColsPerDeg = DEM_NUM_COLS / (abs( DEM_LON_BEG ) + abs( DEM_LON_END ));

	MIN_LAT = double( DEM_LAT_BEG );
	MIN_LON = double( DEM_LON_BEG );

	MAX_LAT = double( DEM_LAT_END ) + (double( 1 ) / double( RowsPerDeg ));
	MAX_LON = double( DEM_LON_END ) - (double( 1 ) / double( ColsPerDeg ));

	ABS_MIN_LAT = fabs( MIN_LAT );
	ABS_MIN_LON = fabs( MIN_LON );

	ABS_MAX_LAT = fabs( MAX_LAT );
	ABS_MAX_LON = fabs( MAX_LON );

	LAT_DENOMINATOR = ABS_MIN_LAT + ABS_MAX_LAT;
	LON_DENOMINATOR = ABS_MIN_LON + ABS_MAX_LON;

	LAT_FACTOR = double( MAX_ROW_INDEX ) / LAT_DENOMINATOR;
	LON_FACTOR = double( MAX_COL_INDEX ) / LON_DENOMINATOR;

	LOCI_BITMASK = 49152; // two most significant bits represent LOCI: 2^15 (16384) + 2^14 (32768) = 49152

	m_fpIn = NULL;
	m_FileArray = NULL;

	m_ReadIntoMemory = FALSE;
}

///////////////////////////////////////////////////////////////////////////////

C_DEM::~C_DEM()
{
	if( m_ReadIntoMemory )
	{
		if( m_FileArray != NULL )
		{
			delete m_FileArray;
		}
	}

	if( m_fpIn != NULL )
	{
		fclose( m_fpIn );
	}
}

///////////////////////////////////////////////////////////////////////////////

DWORD C_DEM::Init( const char *Filename, BOOL ReadIntoMemory )
{
	m_ReadIntoMemory = ReadIntoMemory;

	CFileStatus FileStatus;

	if( !CFile::GetStatus( Filename, FileStatus ) )
	{
		return( GetLastError() );
	}

	ULONGLONG ExpectedFileSize = ULONGLONG( DEM_NUM_ROWS ) * ULONGLONG( DEM_NUM_COLS ) * ULONGLONG( DEM_BYTES_PER_POINT );
	ULONGLONG FileSize = FileStatus.m_size;

	if( ExpectedFileSize != FileSize )
	{
		return CIRAERR_UNEXPECTED_FILESIZE;
	}

	m_fpIn = fopen( Filename, "rb" );
	if( m_fpIn == NULL )
	{
		return( GetLastError() );
	}

	if( ReadIntoMemory )
	{
		ULONGLONG FileArraySize = ULONGLONG( DEM_NUM_ROWS ) * ULONGLONG( DEM_NUM_COLS );

		if( m_FileArray != NULL )
		{
			delete m_FileArray;
		}

		m_FileArray = new unsigned short[FileArraySize];

		// there seems to be a limitation on the number of bytes read at a time using fread (unsigned int max of 4,294,967,295?)
		// since our array is 6,884,352,000 bytes, we will read one half of the array at a time
	 
		unsigned int HalfFileSize = unsigned int( FileSize / ULONGLONG( 2 ) ); // half of the array in bytes is the same as the number of elements

		if( fread( m_FileArray, HalfFileSize, 1, m_fpIn ) != 1 )
		{
			delete m_FileArray;
			m_FileArray = NULL;

			fclose( m_fpIn );
			m_fpIn = NULL;

			return( GetLastError() );
		}

		if( fread( &m_FileArray[FileArraySize/2], HalfFileSize, 1, m_fpIn ) != 1 )
		{
			delete m_FileArray;
			m_FileArray = NULL;

			fclose( m_fpIn );
			m_fpIn = NULL;

			return( GetLastError() );
		}

		fclose( m_fpIn );
		m_fpIn = NULL;
	}

	return( ERROR_SUCCESS );
}

////////////////////////////////////////////////////////////////////////////////

__int64 C_DEM::Round64( double Number )
{
	if( Number < 0.0 )
	{
		Number -= 0.5;
	}
	else
	{
		Number += 0.5;
	}

	return __int64( Number );
}

///////////////////////////////////////////////////////////////////////////////

ULONGLONG C_DEM::Lat2Row( double Lat )
{
	return Round64( (LAT_DENOMINATOR - (Lat + ABS_MAX_LAT)) * LAT_FACTOR );
}

///////////////////////////////////////////////////////////////////////////////

ULONGLONG C_DEM::Lon2Col( double Lon )
{
	ULONGLONG Col = Round64( (Lon + ABS_MIN_LON) * LON_FACTOR );

	if( Col < 0 )
		Col += (MAX_COL_INDEX + 1);
	else if( Col > MAX_COL_INDEX )
		Col -= (MAX_COL_INDEX + 1);

	return Col;
}

///////////////////////////////////////////////////////////////////////////////

double C_DEM::Row2Lat( ULONGLONG Row )
{
	return -(((double( Row ) / LAT_FACTOR) - LAT_DENOMINATOR) + ABS_MAX_LAT);
}

///////////////////////////////////////////////////////////////////////////////

double C_DEM::Col2Lon( ULONGLONG Col )
{
	return (double( Col ) / LON_FACTOR) - ABS_MIN_LON;
}

///////////////////////////////////////////////////////////////////////////////

short C_DEM::GetElevationShort( ULONGLONG Row, ULONGLONG Col, BYTE *LOCI )
{
	if( Row < ULONGLONG( 0 ) ) return DEM_ELEV_ERROR; 
	if( Row > MAX_ROW_INDEX ) return DEM_ELEV_ERROR; 

	if( Col < ULONGLONG( 0 ) ) return DEM_ELEV_ERROR; 
	if( Col > MAX_COL_INDEX ) return DEM_ELEV_ERROR; 

	unsigned short TempVal;

	if( m_ReadIntoMemory )
	{
		ULONGLONG Offset = (Row * COLS_PER_ROW) + Col;
		TempVal = m_FileArray[Offset];
	}
	else
	{
		ULONGLONG Offset = (Row * BYTES_PER_ROW) + (Col * ULONGLONG( DEM_BYTES_PER_POINT ));
		if( _fseeki64( m_fpIn, Offset, SEEK_SET ) != 0 )
		{
			// error
			return DEM_ELEV_ERROR;
		}
		else
		{
			if( fread( &TempVal, DEM_BYTES_PER_POINT, 1, m_fpIn ) != 1 )
			{
				// error
				return DEM_ELEV_ERROR;
			}
		}
	}

	// TempVal is the value read from the file
	// uint16 format
	// 2 most significant bits are the Land Ocean Coast Inland water mask (LOCI).
	// LAND: (when LOCI bits == 00, and hgt > 2000) "DEM_elevation" = Altitude wrt EGS96 GEOID, "Navigation_land_sea_flag" variable set to 1;
	// OCEAN: (when LOCI bits == 00, and hgt < 2000) "DEM_height" = -9999, "Navigation_land_sea_flag" variable set to 2;
	// COASTLINE: DEM_elevation = Altitude wrt EGS96 GEOID, "Navigation_land_sea_flag" variable set to 3;
	// INLAND WATER: DEM_elevation = Altitude wrt EGS96 GEOID, "Navigation_land_sea_flag" variable set to 4;
	// INLAND MIXED: DEM_elevation = Altitude wrt EGS96 GEOID, "Navigation_land_sea_flag" variable set to 5;

	unsigned short TempLOCI = TempVal & LOCI_BITMASK;

	short Elev = short( TempVal - TempLOCI );

	BYTE RawLOCI = BYTE( TempLOCI >> 14 );

	if( Elev > 2000 )
	{
		Elev -= 3000;

		if( LOCI != NULL )
		{
			if( RawLOCI == 0 )
				*LOCI = BYTE( 1 ); // LAND
			else if( RawLOCI == 1 )
				*LOCI = BYTE( 3 ); // COASTLINE
			else if( RawLOCI == 2 )
				*LOCI = BYTE( 4 ); // INLAND WATER
			else if( RawLOCI == 3 )
				*LOCI = BYTE( 5 ); // INLAND MIXED
		}
	}
	else
	{
		Elev = -9999; // OCEAN

		if( LOCI != NULL )
		{
			*LOCI = BYTE( 2 ); // OCEAN
		}
	}

	return Elev;
}

///////////////////////////////////////////////////////////////////////////////

short C_DEM::GetElevationShort( double Lat, double Lon, BYTE *LOCI )
{
	if( m_ReadIntoMemory )
	{
		if( m_FileArray == NULL )
		{
			return DEM_ELEV_ERROR;
		}
	}
	else if( m_fpIn == NULL )
	{
		// error
		return DEM_ELEV_ERROR;
	}

	ULONGLONG Row = Lat2Row( Lat );
	ULONGLONG Col = Lon2Col( Lon );

	return GetElevationShort( Row, Col, LOCI );
}

///////////////////////////////////////////////////////////////////////////////

short C_DEM::GetElevationShort( float Lat, float Lon, BYTE *LOCI )
{
	return GetElevationShort( double( Lat ), double( Lon ), LOCI );
}

///////////////////////////////////////////////////////////////////////////////

float C_DEM::GetElevationFloat( double Lat, double Lon, BYTE *LOCI )
{
	return float( GetElevationShort( Lat, Lon, LOCI ) );
}

///////////////////////////////////////////////////////////////////////////////

float C_DEM::GetElevationFloat( float Lat, float Lon, BYTE *LOCI )
{
	return float( GetElevationShort( double( Lat ), double( Lon ), LOCI ) );
}

