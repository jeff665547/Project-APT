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
//
// ~/apt2-src/util-package/CodeMap.h ---
//
// Revision History:
// Created by David Le on 12/17/2014
// Modified by David Le on 3/18/2014
//

#ifndef __CODEMAP_H__
#define __CODEMAP_H__

#include <map>
#include <vector>
#include <string>
#include <sstream>

enum CallCode {
    CALLCODE_AA = 6,
    CALLCODE_BB = 7,
    CALLCODE_AB = 8,
    CALLCODE_NOCALL = 11,
    CALLCODE_OTV = 12
};

enum CallTranslated {
    TRANSLATED_CALL_AA = 0,
    TRANSLATED_CALL_AB = 1,
    TRANSLATED_CALL_BB = 2,
    TRANSLATED_CALL_NOCALL = -1,
    TRANSLATED_CALL_OTV = -2
};

enum CallFormat {
    FORMAT_NONE = 0,
    FORMAT_CALLCODE = 1,
    FORMAT_BASECALL = 2,
    FORMAT_TRANSLATED = 3,
};

enum CoreColumn {
    CORE_COLUMN_CALLS = 0,
    CORE_COLUMN_CONFIDENCE = 1,
    CORE_COLUMN_LOGRATIO = 2,
    CORE_COLUMN_STRENGTH = 3,
    CORE_COLUMN_A_ALLELE_SIGNAL = 4,
    CORE_COLUMN_B_ALLELE_SIGNAL = 5
};

class CodeMap
{
public:
    static std::map<int, std::string> CallCodes;             // <CallCode, AA:BB:AB:NoCall:OTV>
    static std::map<int, std::string> TranslatedCalls;       // <CallTranslated, AA:BB:AB:NoCall:OTV>
    static std::map<int, CallCode> TranslatedToCallCodes;   // <CallTranslated, 6:7:8:11:12>
    static std::map<int, CallTranslated> CodeToTranslatedCalls; // <CallCode, 0:1:2:-1:-2>
    static std::map<int, std::string> SampleColumns;         // <CoreColumn, call_code:confidence:log_ratio:strength>

public:
    static void init();
    static bool isAlleleCodeValid(std::string& strCode);
    static bool replaceAllelesAB(std::string& strCode, std::string& strAlleleA, std::string& strAlleleB, std::string strSeparator = "");
    static std::string nameCodeToFSBC(std::string strCode, std::string strAlleleA, std::string strAlleleB);
    static std::string nameCodeToVcfByAlt(std::string strCallCode, std::string strAlleleA, std::string strAlleleB, std::string strAltAllele);
    static std::string nameCodeToVcfFSBC(std::string strCode, std::string strAlleleA, std::string strAlleleB);
    static std::string nameCodeToPedFSBC(std::string strCode, std::string strAlleleA, std::string strAlleleB);
    static std::string callCodeToFSBC(CallCode code, std::string strAlleleA, std::string strAlleleB);
    static std::string callCodeToVcf(CallCode code, std::string strAlleleA, std::string strAlleleB, std::string strRefAllele, std::string strAltAllele);
    static std::string callCodeToVcfFSBC(CallCode code, std::string strAlleleA, std::string strAlleleB);
    static std::string callCodeToPedFSBC(CallCode code, std::string strAlleleA, std::string strAlleleB);
    static std::string translatedCallToFSBC(CallTranslated call, std::string strAlleleA, std::string strAlleleB, bool bPedigree = false);
    static std::string codeToTranslatedCall(CallCode code);
    static CallFormat nameFormatToCode(std::string strCode);
    static std::string callFormatToName(CallFormat format);
};

#endif
