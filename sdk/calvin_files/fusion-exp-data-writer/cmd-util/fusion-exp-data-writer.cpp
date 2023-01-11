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

#include <calvin_files/fusion-exp-data-writer/src/ExpressionSignalDetectionWriter.h>
#include <calvin_files/utils/src/StringUtils.h>
#include <calvin_files/converters/utils/src/CmdLine.h>
#include <util/Convert.h>
#include <file/TsvFile/TsvFile.h>

int main(int argc, char* argv[])
{
	try
	{
		CCmdLine cmdLine;
		if (cmdLine.SplitLine(argc, argv) < 1)
			exit(-1);

		std::string dataFileName = cmdLine.GetSafeArgument("-i", 0, "");
		std::string celFileName = cmdLine.GetArgument("-cel", 0);
		std::string chpFileName = celFileName.substr(0, celFileName.length() - 3) + "CHP";
		std::wstring algname = L"SCAN-UPC";
		std::wstring algversion = L"1.0";
		std::wstring progname = L"EC";
		std::wstring progversion = L"1.4";
		std::wstring company = L"Affymetrix";
		std::wstring arrayType = affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(cmdLine.GetArgument("-a", 0));
		ParameterNameValueTypeList algParams;
		ParameterNameValueTypeList summaryParams;
		std::list<std::string> psnames;
		std::list<float> signals;
		std::list<float> pvalues;

		if (dataFileName.length() == 0)
		{
			for (int i=0; i<100; i++)
			{
				psnames.push_back(ToStr(i));
				signals.push_back(i*10.0f);
				pvalues.push_back(i/10.0f);
			}
		}
		else
		{
			std::string psname;
			float sig;
			float pval;
			affx::TsvFile tsv;
			tsv.bind(0, "probeset_id", &psname, affx::TSV_BIND_REQUIRED);
			tsv.bind(0, "signal", &sig, affx::TSV_BIND_REQUIRED);
			tsv.bind(0, "p-value", &pval, affx::TSV_BIND_REQUIRED);
			if (tsv.open(dataFileName) != affx::TSV_OK)
				Err::errAbort("Can't open the input file");

			tsv.headersBegin();
			std::string key;
			std::string val;
			ParameterNameValueType param;
			while (tsv.headersNext(key, val) == affx::TSV_OK)
			{
				param.SetName(affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(key));
				param.SetValueText(affymetrix_calvin_utilities::StringUtils::ConvertMBSToWCS(val));
				ExpressionSignalDetectionWriter::PrependAlgorithmParameterNameTag(param);
				algParams.push_back(param);
			}
			while (tsv.nextLine() == affx::TSV_OK)
			{
				psnames.push_back(psname);
				signals.push_back(sig);
				pvalues.push_back(pval);
			}
			tsv.close();
		}
		
		ExpressionSignalDetectionWriter writer;
		writer.WhiteExpressionSignalDetectionCHPFile(celFileName, chpFileName, algname, algversion, progname, progversion, company, arrayType, algParams, summaryParams, psnames, signals, pvalues);

	}
	catch (std::exception& ex)
	{
		std::cout << ex.what() << std::endl;
	}
	return 0;
}

