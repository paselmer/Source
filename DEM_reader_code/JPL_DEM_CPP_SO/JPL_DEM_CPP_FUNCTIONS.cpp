/* JPL_DEM_CPP_FUNCTIONS.cpp                         */
/* Written by Patrick Selmer. Last modified: 3/15/18 */

/* README

   This particular code was written for Unix compatibily. The Windows-compatible version of this source,
   JPL_DEM_CPP_DLL.cpp does exactly the same task, with any code differences due only to compatibily.
   
   IMPORTANT NOTES:
   
   HOW DO I COMPILE?
   
   g++ -c -Wall -fpic JPL_DEM_CPP_FUNCTIONS.cpp
   g++ -shared -o JPL_DEM_CPP_FUNCTIONS.so JPL_DEM_CPP_FUNCTIONS.o
   
   CODE...
   
   The "extern "C"" needs to be prepended to the "get_DEM_point" function. This ensures
   its name isn't mangled by the compiler. In the Windows version you would also need "__declspec(dllexport)"
   like this ... "extern "C" __declspec(dllexport)." A mangled name btw, means that the calling routine
   won't be able to find the function. This is only an issue for C++, not C. See
   "https://stackoverflow.com/questions/1041866/what-is-the-effect-of-extern-c-in-c" for more info.
   
   There are different header files required for this code versus the Windows version. You can view the headers
   just below. The DOS2UNIX_typedef.h was created to copying to code from the Windows-compatible version.
   The C_DEM_Unix.h file is ALMOST identical to C_DEM.h in Windows.
   
*/   

//#include <afx.h>  ---> I don't have this header. The free version of Visual Studio doesn't come with it. 2/9/18
//#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
//#include <atlstr.h> 
#include "/cpl3/CAMAL/Config/DEM/C_DEM_Unix.h"


///////////////////////////////////* START CODE WRITTEN BY FOLKS AT NASA JPL */
///////////////////////////////////////////////////////////////////////////////
// NOTE: Some of this code was edited by Patrick Selmer (GSFC) for increased
//       compatibility
//
//	C_DEM.cpp
//
//	This class is used to access a digital elevation map (DEM) in order to
//  get the elevation at a specified latitude and longitude
//
//	Notes:
//		Returned elevation of 9999 indicates an error
//

///////////////////////////////////////////////////////////////////////////////

C_DEM::C_DEM()
{
	COLS_PER_ROW = ULONGLONG(DEM_NUM_COLS);
	BYTES_PER_ROW = ULONGLONG(DEM_NUM_COLS) * ULONGLONG(DEM_BYTES_PER_POINT);

	MAX_ROW_INDEX = ULONGLONG(DEM_NUM_ROWS - 1);
	MAX_COL_INDEX = ULONGLONG(DEM_NUM_COLS - 1);

	int RowsPerDeg = DEM_NUM_ROWS / (abs(DEM_LAT_BEG) + abs(DEM_LAT_END));
	int ColsPerDeg = DEM_NUM_COLS / (abs(DEM_LON_BEG) + abs(DEM_LON_END));

	MIN_LAT = double(DEM_LAT_BEG);
	MIN_LON = double(DEM_LON_BEG);

	MAX_LAT = double(DEM_LAT_END) + (double(1) / double(RowsPerDeg));
	MAX_LON = double(DEM_LON_END) - (double(1) / double(ColsPerDeg));

	ABS_MIN_LAT = fabs(MIN_LAT);
	ABS_MIN_LON = fabs(MIN_LON);

	ABS_MAX_LAT = fabs(MAX_LAT);
	ABS_MAX_LON = fabs(MAX_LON);

	LAT_DENOMINATOR = ABS_MIN_LAT + ABS_MAX_LAT;
	LON_DENOMINATOR = ABS_MIN_LON + ABS_MAX_LON;

	LAT_FACTOR = double(MAX_ROW_INDEX) / LAT_DENOMINATOR;
	LON_FACTOR = double(MAX_COL_INDEX) / LON_DENOMINATOR;

	LOCI_BITMASK = 49152; // two most significant bits represent LOCI: 2^15 (16384) + 2^14 (32768) = 49152

	m_fpIn = NULL;
	m_FileArray = NULL;

	m_ReadIntoMemory = FALSE;
}

///////////////////////////////////////////////////////////////////////////////

C_DEM::~C_DEM()
{
	if (m_ReadIntoMemory)
	{
		if (m_FileArray != NULL)
		{
			delete m_FileArray;
		}
	}

	if (m_fpIn != NULL)
	{
		fclose(m_fpIn);
	}
}

///////////////////////////////////////////////////////////////////////////////

DWORD C_DEM::Init(const char *Filename, BOOL ReadIntoMemory)
{
	m_ReadIntoMemory = ReadIntoMemory;


	/* This error check removed 2/9/18 because I don't want to pay for the
	version of Visual Studio that uses it
	CFileStatus FileStatus;

	if (!CFile::GetStatus(Filename, FileStatus))
	{
		return(DEM_FILE_READ_ERROR);
	}
	*/

	ULONGLONG ExpectedFileSize = ULONGLONG(DEM_NUM_ROWS) * ULONGLONG(DEM_NUM_COLS) * ULONGLONG(DEM_BYTES_PER_POINT);
	/*ULONGLONG FileSize = FileStatus.m_size;   removed 2/9/18

	if (ExpectedFileSize != FileSize)
	{
		return CIRAERR_UNEXPECTED_FILESIZE;
	}
	*/

	m_fpIn = fopen(Filename, "rb");
	if (m_fpIn == NULL)
	{
		return(DEM_FILE_READ_ERROR);
	}

	if (ReadIntoMemory)
	{
		ULONGLONG FileArraySize = ULONGLONG(DEM_NUM_ROWS) * ULONGLONG(DEM_NUM_COLS);

		if (m_FileArray != NULL)
		{
			delete m_FileArray;
		}

		m_FileArray = new unsigned short[FileArraySize];

		// there seems to be a limitation on the number of bytes read at a time using fread (unsigned int max of 4,294,967,295?)
		// since our array is 6,884,352,000 bytes, we will read one half of the array at a time

		// The following single line was modified to use ExpectedFileSize instead of FileSize because I didn't have the right header 2/9/18
		unsigned int HalfFileSize = UINT(ExpectedFileSize / ULONGLONG(2)); // half of the array in bytes is the same as the number of elements

		if (fread(m_FileArray, HalfFileSize, 1, m_fpIn) != 1)
		{
			delete m_FileArray;
			m_FileArray = NULL;

			fclose(m_fpIn);
			m_fpIn = NULL;

			return(DEM_FILE_READ_ERROR);
		}

		if (fread(&m_FileArray[FileArraySize / 2], HalfFileSize, 1, m_fpIn) != 1)
		{
			delete m_FileArray;
			m_FileArray = NULL;

			fclose(m_fpIn);
			m_fpIn = NULL;

			return(DEM_FILE_READ_ERROR);
		}

		fclose(m_fpIn);
		m_fpIn = NULL;
	}

	return(ERROR_SUCCESS);
}

////////////////////////////////////////////////////////////////////////////////

__int64 C_DEM::Round64(double Number)
{
	if (Number < 0.0)
	{
		Number -= 0.5;
	}
	else
	{
		Number += 0.5;
	}

	return __int64(Number);
}

///////////////////////////////////////////////////////////////////////////////

ULONGLONG C_DEM::Lat2Row(double Lat)
{
	return Round64((LAT_DENOMINATOR - (Lat + ABS_MAX_LAT)) * LAT_FACTOR);
}

///////////////////////////////////////////////////////////////////////////////

ULONGLONG C_DEM::Lon2Col(double Lon)
{
	ULONGLONG Col = Round64((Lon + ABS_MIN_LON) * LON_FACTOR);

	if (Col < 0)
		Col += (MAX_COL_INDEX + 1);
	else if (Col > MAX_COL_INDEX)
		Col -= (MAX_COL_INDEX + 1);

	return Col;
}

///////////////////////////////////////////////////////////////////////////////

double C_DEM::Row2Lat(ULONGLONG Row)
{
	return -(((double(Row) / LAT_FACTOR) - LAT_DENOMINATOR) + ABS_MAX_LAT);
}

///////////////////////////////////////////////////////////////////////////////

double C_DEM::Col2Lon(ULONGLONG Col)
{
	return (double(Col) / LON_FACTOR) - ABS_MIN_LON;
}

///////////////////////////////////////////////////////////////////////////////

short C_DEM::GetElevationShort(ULONGLONG Row, ULONGLONG Col, BYTE *LOCI)
{
	if (Row < ULONGLONG(0)) return DEM_ELEV_ERROR;
	if (Row > MAX_ROW_INDEX) return DEM_ELEV_ERROR;

	if (Col < ULONGLONG(0)) return DEM_ELEV_ERROR;
	if (Col > MAX_COL_INDEX) return DEM_ELEV_ERROR;

	unsigned short TempVal;

	if (m_ReadIntoMemory)
	{
		ULONGLONG Offset = (Row * COLS_PER_ROW) + Col;
		TempVal = m_FileArray[Offset];
	}
	else
	{
		ULONGLONG Offset = (Row * BYTES_PER_ROW) + (Col * ULONGLONG(DEM_BYTES_PER_POINT));
		if (fseeko(m_fpIn, Offset, SEEK_SET) != 0)
		{
			// error
			return DEM_ELEV_ERROR;
		}
		else
		{
			if (fread(&TempVal, DEM_BYTES_PER_POINT, 1, m_fpIn) != 1)
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

	short Elev = short(TempVal - TempLOCI);

	BYTE RawLOCI = BYTE(TempLOCI >> 14);

	if (Elev > 2000)
	{
		Elev -= 3000;

		if (LOCI != NULL)
		{
			if (RawLOCI == 0)
				*LOCI = BYTE(1); // LAND
			else if (RawLOCI == 1)
				*LOCI = BYTE(3); // COASTLINE
			else if (RawLOCI == 2)
				*LOCI = BYTE(4); // INLAND WATER
			else if (RawLOCI == 3)
				*LOCI = BYTE(5); // INLAND MIXED
		}
	}
	else
	{
		Elev = -9999; // OCEAN

		if (LOCI != NULL)
		{
			*LOCI = BYTE(2); // OCEAN
		}
	}

	return Elev;
}

///////////////////////////////////////////////////////////////////////////////

short C_DEM::GetElevationShort(double Lat, double Lon, BYTE *LOCI)
{
	if (m_ReadIntoMemory)
	{
		if (m_FileArray == NULL)
		{
			return DEM_ELEV_ERROR;
		}
	}
	else if (m_fpIn == NULL)
	{
		// error
		return DEM_ELEV_ERROR;
	}

	ULONGLONG Row = Lat2Row(Lat);
	ULONGLONG Col = Lon2Col(Lon);

	return GetElevationShort(Row, Col, LOCI);
}

///////////////////////////////////////////////////////////////////////////////

short C_DEM::GetElevationShort(float Lat, float Lon, BYTE *LOCI)
{
	return GetElevationShort(double(Lat), double(Lon), LOCI);
}

///////////////////////////////////////////////////////////////////////////////

float C_DEM::GetElevationFloat(double Lat, double Lon, BYTE *LOCI)
{
	return float(GetElevationShort(Lat, Lon, LOCI));
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////* END CODE WRITTEN BY FOLKS AT NASA JPL */



extern "C" DEM_outputs get_DEM_point(const char *DEM_file_name, double inlat, double inlon)
{

	/* This function returns a structure composed of the DEM elevation in meters and the 
	   surface type. This function accepts a null-terminated C-style character string file name of the DEM file
	   as well as the lat and lon in degrees.

	   Written by Patrick Selmer
	*/

	C_DEM DEM;
	DEM_outputs output;
	output.Elev_rtn = -99; // Initialize to -99/99
	output.LOCI_rtn = 99;

	unsigned long Status = DEM.Init(DEM_file_name, FALSE);
	if (Status == ERROR_SUCCESS)
	{
		unsigned char LOCI;

		//float Lat = float(40.5592);
		//float Lon = float(-105.0781);
		double Lat = inlat;
		double Lon = inlon;

		float Elev = DEM.GetElevationFloat(Lat, Lon, &LOCI);

		// NOTE: the raw LOCI flag extracted directly from the DEM file is different than the LOCI flag that is returned from the read routine.
		// the raw LOCI flag is a system flag designed to save space (see DEM_description.txt for more info) and the LOCI flag returned here is the finalized science version of the flag.
		
		/* All I need is LOCI, not a "CString" type description... [3/14/18]

		CString LOCI_Description;

		if (LOCI == 1)
			LOCI_Description = CString("land");
		else if (LOCI == 2)
			LOCI_Description = CString("ocean");
		else if (LOCI == 3)
			LOCI_Description = CString("coast");
		else if (LOCI == 4)
			LOCI_Description = CString("inland waters");
		else if (LOCI == 5)
			LOCI_Description = CString("inland mixed");

		//printf("Elev = %f | LOCI = %d (%s)\n", Elev, LOCI, LOCI_Description);
		*/

		output.Elev_rtn = Elev;
		output.LOCI_rtn = LOCI;
		return output;
	}
	else {
		printf("C++ code failed for some reason.\n");
		printf("Error Status code is: %ld\n", Status);
		return output;
	}

}
