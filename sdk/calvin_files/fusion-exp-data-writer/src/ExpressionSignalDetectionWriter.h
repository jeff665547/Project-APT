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

#pragma once

#include <calvin_files/writers/src/CalvinCHPQuantificationDetectionFileWriter.h>
#include <calvin_files/data/src/ProbeSetQuantificationDetectionData.h>

using namespace affymetrix_calvin_io;
using namespace affymetrix_calvin_data;

/*! This class provides an easy to use wrapper to creating expression signal/detection CHP files. */
class ExpressionSignalDetectionWriter
{
private:
	/*! The data object. */
	CHPQuantificationDetectionData *pdata;

public:
	/*! Constructor */
	ExpressionSignalDetectionWriter();

	/*! Destructor */
	~ExpressionSignalDetectionWriter();

	/*! Prepends the required tag to the name of an algorithm parameter.
	 * @param param The parameter to modify.
	 */
	static void PrependAlgorithmParameterNameTag(ParameterNameValueType &param);
	
	/*! Prepends the required tag to the name of an algorithm option parameter.
	 * @param param The parameter to modify.
	 */
	static void PrependAlgorithmParameterOptionNameTag(ParameterNameValueType &param);

private:
	/*! Add the CEL header to the CHP object.
	 * @param celFileName The name of the CEL file (full path)
	 */
	void AddCELHeader(const std::string &celFileName);
	
	/*! Add the algorithm information to the CHP object.
	 * @param algname The algorithm name.
	 * @param algversion The algorithm version
	 * @param progname The program name
	 * @param progversion The program version
	 * @param progcompany The company name.
	 */
	void AddAlgorithmInformation(const wchar_t *algname, const wchar_t *algversion, const wchar_t *progname, const wchar_t *progversion, const wchar_t *progcompany);
	
	/*! Add the array information to the CHP object.
	 * @param nsets The number of probe sets
	 * @param maxLen The maximum length of the probe set names
	 * @param arrayType The probe array type
	 */
	void AddArrayInfo(int nsets, int maxLen, const wchar_t *arrayType);
	
	/*! Add the parameters to the CHP object.
	 * @param params The parameters to add
	 * @param storeAsParams True if stored as a parameter, false if stored as a summary metric.
	 */
	void AddParameters(const ParameterNameValueTypeList &params, bool storeAsParams);

public:
	/*! Add the algorithm information to the CHP object.
	 * @param celFileName The name of the CEL file (full path)
	 * @param chpFileName The name of the CHP file to write (full path)
	 * @param algname The algorithm name.
	 * @param algversion The algorithm version
	 * @param progname The program name
	 * @param progversion The program version
	 * @param progcompany The company name.
 	 * @param arrayType The probe array type
	 * @param algParams The algorithm parameters
	 * @param summaryParams The summary metrics
	 * @param psnames The probe set names
	 * @param signals The signal values
	 * @param pvalues The p-values
	 */
	void WhiteExpressionSignalDetectionCHPFile(
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
			const std::list<float> &pvalues);
};
