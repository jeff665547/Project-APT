////////////////////////////////////////////////////////////////
//
// Copyright (C) 2016 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////
#include <string>
#include <vector>

#include "AGCC_FileUtilities.h"

using namespace std;

#include "calvin_files/parsers/src/CelFileReader.h"
#include "calvin_files/parsers/src/GenericFileReader.h"
using namespace affymetrix_calvin_io;

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef UINT
#define UINT unsigned int
#endif

#ifndef WORD
#define WORD unsigned short
#endif

#ifndef errno_t
#define errno_t int
#endif

#ifndef byte
#define byte unsigned char
#endif


bool copy_file (const char *old_file, const char *new_file)
{
	FILE *stream = fopen (old_file,"rb");
	if (stream == 0)
		return false;

	FILE *stream_out = fopen (new_file,"wb");
	if (stream_out == 0)
	{
		fclose (stream);
		return false;
	}

	unsigned long sz = 0;
	while (1)
	{
		int c = fgetc (stream);
		if (c == EOF)
			break;

		sz++;

		fputc (c,stream_out);
	}

	fclose (stream);
	fclose (stream_out);

	return true;
}

static DataSet *GetDataSet (wstring name, wstring group, GenericData &hdr_data)
{
	WStringVector groupNames;
	hdr_data.DataGroupNames(groupNames);
	int32_t ngroups = (int32_t)groupNames.size();
	for (int32_t igroup=0; igroup<ngroups; igroup++)
	{
		wstring check = groupNames[igroup];
		if (check != group)
			continue;

		int32_t nsets = hdr_data.DataSetCnt (igroup);

		for (int32_t iset=0; iset<nsets; iset++)
		{
			DataSet *data_set = hdr_data.DataSet (igroup,iset);

			wstring data_set_name = data_set->Header ().GetName ();
			if (data_set_name == name)
			{
				return data_set;
			}
		}
	}

	return 0;
}

static int Fusion_GetDataSetOffsetAndSize
(
	const char *file_name,
	const char *wavelength,

	int &offset,
	int &size,

	char *err_msg,
	int err_msg_alloc_sz
	)
{
	GenericFileReader reader;

	string fname = file_name;
	reader.SetFilename (fname);

	GenericData hdr_data;
	reader.ReadHeader (hdr_data);

	GenericDataHeader *hdr = hdr_data.Header ().GetGenericDataHdr ();

	ParameterNameValueType cel_cols,cel_rows;	
	if (hdr->FindNameValParam(L"affymetrix-cel-cols",cel_cols) == false || 
		hdr->FindNameValParam(L"affymetrix-cel-cols",cel_rows) == false)
	{
		strncpy  (err_msg,"Invalid cel file",sizeof(err_msg_alloc_sz)-1);
		return FALSE;
	}

	string err_s;

	while (1)
	{
		size_t len = strlen (wavelength)+1;
		wchar_t *_wavelength = new wchar_t[len];
		for (UINT i=0; i<len; i++)
		{
			_wavelength[i] = (WORD)wavelength[i];
		}

		wstring group_name = _wavelength;
		delete [] _wavelength;

		DataGroupHeader *data_group_header = hdr_data.Header ().FindDataGroupHeader (group_name);

		if (data_group_header == 0)
		{
			err_s = "Missing data group ";
			err_s += wavelength;
			break;
		}

		DataSet *data_set = GetDataSet (L"Intensity", group_name, hdr_data);

		if (data_set == 0)
		{
			err_s = "Missing data set \"Intensity\"";
			break;
		}

		offset = data_set->Header ().GetDataStartFilePos ();
		size = data_set->Header ().GetDataSize ();

		break;
	}

	if (strlen (err_s.c_str ()))
	{
		strncpy (err_msg,err_s.c_str (),sizeof(err_msg_alloc_sz)-1);
		return FALSE;
	}

	return TRUE;
}

int Fusion_UpdateCelFileIntensityValues
(
	const char *file_name,
	const char *wavelength,


	std::vector<float>& intensity_array,
	int array_sz,

	char *err_msg,
	int err_msg_alloc_sz
)
{
	int offset=0,size=0;
	if (Fusion_GetDataSetOffsetAndSize (file_name,wavelength,offset,size,err_msg,err_msg_alloc_sz) == FALSE)
		return FALSE;

	string err_s;

	while (1)
	{
		FILE *stream = fopen (file_name,"r+b");
		if (!stream)
		{
			err_s = "Can't open file: ";
			err_s += file_name;
			err_s += " for updating.";
			break;
		}

		if (fseek (stream,offset,SEEK_SET))
		{
			err_s = "Problem seeking within file: ";
			err_s += file_name;
			break;
		}

		int array_sz_bytes = array_sz * sizeof(float);
		if (array_sz_bytes != size)
		{
			err_s = "Given array size does not match data set size in file: ";
			err_s += file_name;
			fclose (stream);
			break;
		}

		int i=0;
		for (i=0; i<array_sz; i++)
		{
			union 
			{
				float f;
				byte b[4];
			} source,target;

			source.f = intensity_array[i];
			target.b[0] = source.b[3];
			target.b[1] = source.b[2];
			target.b[2] = source.b[1];
			target.b[3] = source.b[0];
			
			if (fwrite (&target.b[0],4,1,stream) != 1)
				break;
		}

		if (i != array_sz)
		{
			err_s = "Problem writing data to file: ";
			err_s += file_name;
			fclose (stream);
			break;
		}

		fclose (stream);

		break;
	}

	if (strlen (err_s.c_str ()))
	{
		strncpy (err_msg,err_s.c_str (),sizeof(err_msg_alloc_sz)-1);
		return FALSE;
	}

	return TRUE;
}
