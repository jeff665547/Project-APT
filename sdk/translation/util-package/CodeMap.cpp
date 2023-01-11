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
// ~/apt2-src/util-package/CodeMap.cpp ---
//
// Revision History:
// Created by David Le on 12/17/2014
// Modified by David Le on 3/18/2014
//

#include "CodeMap.h"

#include "util/Convert.h"

using namespace std;

map<int, string> CodeMap::CallCodes;
map<int, string> CodeMap::TranslatedCalls;
map<int, CallCode> CodeMap::TranslatedToCallCodes;
map<int, CallTranslated> CodeMap::CodeToTranslatedCalls;
map<int, string> CodeMap::SampleColumns;

void CodeMap::init()
{
    CodeMap::CallCodes[CALLCODE_AA] = "AA";
    CodeMap::CallCodes[CALLCODE_BB] = "BB";
    CodeMap::CallCodes[CALLCODE_AB] = "AB";
    CodeMap::CallCodes[CALLCODE_NOCALL] = "NoCall";
    CodeMap::CallCodes[CALLCODE_OTV] = "OTV";

    CodeMap::TranslatedCalls[TRANSLATED_CALL_AA] = "AA";
    CodeMap::TranslatedCalls[TRANSLATED_CALL_AB] = "AB";
    CodeMap::TranslatedCalls[TRANSLATED_CALL_BB] = "BB";
    CodeMap::TranslatedCalls[TRANSLATED_CALL_NOCALL] = "NoCall";
    CodeMap::TranslatedCalls[TRANSLATED_CALL_OTV] = "OTV";

    CodeMap::TranslatedToCallCodes[TRANSLATED_CALL_AA] = CALLCODE_AA;
    CodeMap::TranslatedToCallCodes[TRANSLATED_CALL_AB] = CALLCODE_AB;
    CodeMap::TranslatedToCallCodes[TRANSLATED_CALL_BB] = CALLCODE_BB;
    CodeMap::TranslatedToCallCodes[TRANSLATED_CALL_NOCALL] = CALLCODE_NOCALL;
    CodeMap::TranslatedToCallCodes[TRANSLATED_CALL_OTV] = CALLCODE_OTV;

    CodeMap::CodeToTranslatedCalls[CALLCODE_AA] = TRANSLATED_CALL_AA;
    CodeMap::CodeToTranslatedCalls[CALLCODE_AB] = TRANSLATED_CALL_AB;
    CodeMap::CodeToTranslatedCalls[CALLCODE_BB] = TRANSLATED_CALL_BB;
    CodeMap::CodeToTranslatedCalls[CALLCODE_NOCALL] = TRANSLATED_CALL_NOCALL;
    CodeMap::CodeToTranslatedCalls[CALLCODE_OTV] = TRANSLATED_CALL_OTV;

    CodeMap::SampleColumns[CORE_COLUMN_CALLS] = "call_code";
    CodeMap::SampleColumns[CORE_COLUMN_CONFIDENCE] = "confidence";
    CodeMap::SampleColumns[CORE_COLUMN_LOGRATIO] = "log_ratio";
    CodeMap::SampleColumns[CORE_COLUMN_STRENGTH] = "strength";
    CodeMap::SampleColumns[CORE_COLUMN_A_ALLELE_SIGNAL] = "a_allele_signal";
    CodeMap::SampleColumns[CORE_COLUMN_B_ALLELE_SIGNAL] = "b_allele_signal";
}

bool CodeMap::isAlleleCodeValid(string& strCode)
{
    return (strCode == "AA" || strCode == "BB" || strCode == "AB");
}

bool CodeMap::replaceAllelesAB(string& strCode, string& strAlleleA, string& strAlleleB, string strSeparator)
{
    if (strAlleleA.length() >= 1 && strAlleleB.length() >= 1) {
        string strReplacedCode = "";
        for (int i = 0; i < strCode.length(); i++) {
            if (strCode[i] == 'A') {
                strReplacedCode += strAlleleA;
            }
            else if (strCode[i] == 'B') {
                strReplacedCode += strAlleleB;
            }
            else {
                strReplacedCode += strCode[i];
            }
            if (i < strCode.length() - 1) {
                strReplacedCode += strSeparator;
            }
        }
        strCode = strReplacedCode;
        return true;
    }
    else {
        return false;
    }
}

string CodeMap::nameCodeToFSBC(string strCode, string strAlleleA, string strAlleleB)
{
    // FSBC: forward strand base call for TXT format
    string strPrevCode = strCode;
    if (CodeMap::isAlleleCodeValid(strCode)) {
        if (CodeMap::replaceAllelesAB(strCode, strAlleleA, strAlleleB)) {
            return strCode;
        }
        else {
            return "00";
        }
    }
    else {
        return strCode;
    }
}

string CodeMap::nameCodeToVcfByAlt(string strCallCode, string strAlleleA, string strAlleleB, string strAltAllele)
{
    vector<string> tokens;
    istringstream line(strAltAllele);
    string token;
    while (getline(line, token, '/')) {
        tokens.push_back(token);
    }

    if ((int)tokens.size() < 2) {
        return "./.";
    }

    string vcfAlleleA = ".";
    string vcfAlleleB = ".";
    bool refAlleleFound = false;
    bool altAlleleFound = false;
    for (int i = 0; i < (int)tokens.size(); i++) {
        if (tokens[i] == strAlleleA) {
            if (!refAlleleFound) {
                vcfAlleleA = Convert::toString(i + 1);
            }
            else if (!altAlleleFound) {
                vcfAlleleB = Convert::toString(i + 1);
            }
        }
        else if (tokens[i] == strAlleleB) {
            if (!refAlleleFound) {
                vcfAlleleB = Convert::toString(i + 1);
            }
            else if (!altAlleleFound) {
                vcfAlleleA = Convert::toString(i + 1);
            }
        }
    }

    if (vcfAlleleA == "." || vcfAlleleB == ".") {
        return "./.";
    }
    else {
        return CodeMap::nameCodeToVcfFSBC(strCallCode, vcfAlleleA, vcfAlleleB);
    }
}

string CodeMap::nameCodeToVcfFSBC(string strCode, string strAlleleA, string strAlleleB)
{
    // FSBC: forward strand base call for plink VCF format
    string strPrevCode = strCode;
    if (isAlleleCodeValid(strCode)) {
        if (CodeMap::replaceAllelesAB(strCode, strAlleleA, strAlleleB, "/")) {
            if (strCode[0] > strCode[2]) {
                char chTemp = strCode[0];
                strCode[0] = strCode[2];
                strCode[2] = chTemp;
            }
            return strCode;
        }
        else {
            return "./.";
        }
    }
    else {
        return "./.";
    }
}

string CodeMap::nameCodeToPedFSBC(string strCode, string strAlleleA, string strAlleleB)
{
    // FSBC: forward strand base call for plink PED format
    string strPrevCode = strCode;
    if (isAlleleCodeValid(strCode)) {
        if (CodeMap::replaceAllelesAB(strCode, strAlleleA, strAlleleB, " ")) {
            return strCode;
        }
        else {
            return "0 0";
        }
    }
    else {
        return "0 0";
    }
}

string CodeMap::callCodeToFSBC(CallCode code, string strAlleleA, string strAlleleB)
{
    string strCallCode = CodeMap::CallCodes[code];
    string strFSBC = "---"; // NoCall or OTV
    if (CodeMap::isAlleleCodeValid(strCallCode)) {
        strFSBC = CodeMap::nameCodeToFSBC(strCallCode, strAlleleA, strAlleleB);
        if (strFSBC == "" || strFSBC == "00") {
            strFSBC = "---"; // true if one of the alleles is empty (similar to GTC)
        }
    }
    return strFSBC;
}

string CodeMap::callCodeToVcf(CallCode code, string strAlleleA, string strAlleleB, string strRefAllele, string strAltAllele)
{
    string strCallCode = CodeMap::CallCodes[code];
    if (strRefAllele == "" || strAltAllele == "") {
        return "./.";
    }
    else if (strRefAllele == ".") {
        return CodeMap::nameCodeToVcfByAlt(strCallCode, strAlleleA, strAlleleB, strAltAllele);
    }
    else if (strRefAllele == strAlleleA) {
        return CodeMap::nameCodeToVcfFSBC(strCallCode, "0", "1");
    }
    else if (strRefAllele == strAlleleB) {
        return CodeMap::nameCodeToVcfFSBC(strCallCode, "1", "0");
    }
    else {
        return "./.";
    }
}

string CodeMap::callCodeToVcfFSBC(CallCode code, string strAlleleA, string strAlleleB)
{
    string strCallCode = CodeMap::CallCodes[code];
    return CodeMap::nameCodeToVcfFSBC(strCallCode, strAlleleA, strAlleleB);
}

string CodeMap::callCodeToPedFSBC(CallCode code, string strAlleleA, string strAlleleB)
{
    string strCallCode = CodeMap::CallCodes[code];
    return CodeMap::nameCodeToPedFSBC(strCallCode, strAlleleA, strAlleleB);
}

string CodeMap::translatedCallToFSBC(CallTranslated call, string strAlleleA, string strAlleleB, bool bPedigree)
{
    string strCall = CodeMap::TranslatedCalls[call];
    string strFSBC;
    if (bPedigree) {
        strFSBC = CodeMap::nameCodeToPedFSBC(strCall, strAlleleA, strAlleleB);
    }
    else {
        strFSBC = CodeMap::nameCodeToFSBC(strCall, strAlleleA, strAlleleB);
    }
    if (bPedigree && strFSBC == strCall) {
        return "0 0";
    }
    return strFSBC;
}

string CodeMap::codeToTranslatedCall(CallCode code)
{
    return Convert::toString(CodeMap::CodeToTranslatedCalls[code]);
}

CallFormat CodeMap::nameFormatToCode(string strCode)
{
    if (strCode == "call_code") { return FORMAT_CALLCODE; }
    if (strCode == "base_call") { return FORMAT_BASECALL; }
    if (strCode == "translated") { return FORMAT_TRANSLATED; }
    return FORMAT_NONE;
};

string CodeMap::callFormatToName(CallFormat format)
{
    if (format == FORMAT_CALLCODE) { return "call_code"; }
    if (format == FORMAT_BASECALL) { return "base_call"; }
    if (format == FORMAT_TRANSLATED) { return "translated"; }
    return "undefined";
};
