////////////////////////////////////////////////////////////////
//
// Copyright (C) 2016 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// (C) Copyright 2014, Affymetrix, Inc.
// All rights reserved. Confidential. Except as pursuant
// to a written agreement with Affymetrix, this software may
// not be used or distributed. This software may be covered
// by one or more patents.
//
// "GeneChip", "Affymetrix" and the Affymetrix logos, and
// Affymetrix user interfaces are trademarks used by Affymetrix.
/////////////////////////////////////////////////////////////////////

#include "ExpressionSignalDetectionWriter.h"
#include <calvin_files/fusion/src/FusionCELData.h>
#include <calvin_files/data/src/GenericData.h>
#include <util/Err.h>

using namespace affymetrix_fusion_io;

ExpressionSignalDetectionWriter::ExpressionSignalDetectionWriter()
{
	pdata = NULL;
}

ExpressionSignalDetectionWriter::~ExpressionSignalDetectionWriter()
{
	pdata = NULL;
}

static void PrepentTag(ParameterNameValueType &param, const wchar_t *tag)
{
	if (param.GetName().find(tag) == -1)
		param.SetName(tag + param.GetName());
}

void ExpressionSignalDetectionWriter::PrependAlgorithmParameterNameTag(ParameterNameValueType &param)
{
	const wchar_t* tag = L"apt-";
	PrepentTag(param, tag);
}

void ExpressionSignalDetectionWriter::PrependAlgorithmParameterOptionNameTag(ParameterNameValueType &param)
{
	const wchar_t* tag = L"apt-opt-";
	PrepentTag(param, tag);
}

void ExpressionSignalDetectionWriter::AddCELHeader(const std::string &celFileName)
{
	FusionCELData cel;
    cel.SetFileName(celFileName.c_str());
    if (!cel.ReadHeader())
		Err::errAbort("Unable to read CEL file: " + celFileName);
	GenericData *gdata = cel.GetGenericData();
	if (gdata != NULL)
		pdata->GetFileHeader()->GetGenericDataHdr()->AddParent(*gdata->Header().GetGenericDataHdr());
	cel.Close();
}

void ExpressionSignalDetectionWriter::AddAlgorithmInformation
(
	const wchar_t *algname,
	const wchar_t *algversion,
	const wchar_t *progname,
	const wchar_t *progversion,
	const wchar_t *progcompany
)
{
	pdata->SetAlgName(algname);
	pdata->SetAlgVersion(algversion);
	
	ParameterNameValueType param;
	const wchar_t *names[] = { L"program-name", L"program-version", L"program-company" };
	const wchar_t *values[] = { progname, progversion, progcompany };
	for (int i=0; i<sizeof(names) / sizeof(wchar_t*); i++)
	{
		param.SetName(names[i]);
		param.SetValueText(values[i]);
		pdata->GetGenericData().Header().GetGenericDataHdr()->AddNameValParam(param);
	}
}

void ExpressionSignalDetectionWriter::AddArrayInfo(int nsets, int maxLen, const wchar_t *arrayType)
{
	pdata->SetEntryCount(nsets, maxLen); 
	pdata->SetArrayType(arrayType);
}

void ExpressionSignalDetectionWriter::AddParameters(const ParameterNameValueTypeList &params, bool storeAsParams)
{
	if (storeAsParams == true)
		pdata->AddAlgParams(params);
	else
		pdata->AddSummaryParams(params);
}

void ExpressionSignalDetectionWriter::WhiteExpressionSignalDetectionCHPFile
(
	const std::string &celFileName,
	const std::string &chpFileName,
	const std::wstring &algname,
	const std::wstring &algversion,
	const std::wstring &progname,
	const std::wstring &progversion,
	const std::wstring &company,
	const std::wstring &arrayType,
	const ParameterNameValueTypeList &algParams,
	const ParameterNameValueTypeList &summaryParams,
	const std::list<std::string> &psnames,
	const std::list<float> &signals,
	const std::list<float> &pvalues
)
{
	int nsets = (int) psnames.size();
	int maxLen = 0;
	for (std::list<std::string>::const_iterator it=psnames.begin(); it!=psnames.end(); it++)
		maxLen = max(maxLen, (int) it->length());

	CHPQuantificationDetectionData data(chpFileName);
	pdata = &data;
	AddCELHeader(celFileName);
	AddAlgorithmInformation(algname.c_str(), algversion.c_str(), progname.c_str(), progversion.c_str(), company.c_str());
	AddArrayInfo(nsets, maxLen, arrayType.c_str());
	AddParameters(algParams, true);
	AddParameters(summaryParams, true);
				
    ProbeSetQuantificationDetectionData entry;
    CHPQuantificationDetectionFileWriter writer(data);
    writer.SeekToDataSet();
	
	int idx=0;
	std::list<std::string>::const_iterator psit = psnames.begin();
	std::list<float>::const_iterator sigit = signals.begin();
	std::list<float>::const_iterator pvalit = pvalues.begin();
    while (psit != psnames.end())
	{
		entry.id = idx++;
		entry.name = *psit++;
        entry.quantification = *sigit++;
		entry.pvalue = *pvalit++;
        writer.WriteEntry(entry);
    }
}
