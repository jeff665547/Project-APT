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
// ~/apt2-src/util-package/CallTable.cpp ---
//
// Revision History:
// Created by David Le on 11/30/2015
//

#include "CallTable.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/Util.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string.h>
#include <memory>

using namespace std;

CallTable::CallTable()
{
    reset();
}

CallTable::~CallTable()
{
}

void CallTable::reset()
{
    m_strErrMsg = "";
    m_strIndexFile = "";

    // default to bi-allelic
    m_iMaxAlleles = 2; 
    m_iMaxCNStates = 2;
    m_iNextRowIndex = 0;
}

bool CallTable::readFile(string strIndexFile)
{
    ifstream file;
    file.open(strIndexFile.c_str());
    if (!file.is_open()) {
        setErrMsg("Failed to open index file for reading.\tFile: " + strIndexFile);
        return false;
    }

    reset();
    char szLine[MAX_LINE];
    while (file.getline(szLine, MAX_LINE)) {
        string strLine = szLine;
        Util::trimString(strLine);
        if (strLine.find("#%%max-alleles=") == 0) {
            string strMaxAlleles = strLine.substr(15);
            Util::trimString(strMaxAlleles);
            m_iMaxAlleles = Convert::toInt(strMaxAlleles);
        }
        else if (strLine.find("#%%max-cn-states=") == 0) {
            string strMaxCNStates = strLine.substr(17);
            Util::trimString(strMaxCNStates);
            m_iMaxCNStates = Convert::toInt(strMaxCNStates);
        }
        else if (strLine.find("#%%call-code-") == 0) {
            vector<string> vTokens, vParts;
            Util::chopString(strLine, "=", vTokens);
            if (vTokens.size() == 2) {
                string strField = vTokens[1];
                Util::trimString(strField);
                Util::chopString(strField, ":", vParts);
                if (vParts.size() == 3) {
                    string strCallName = vParts[0]; Util::trimString(strCallName);
                    string strCallCode = vParts[1]; Util::trimString(strCallCode);
                    string strCNState = vParts[2]; Util::trimString(strCNState);
                    int iCallCode, iCNState;
                    iCallCode = Convert::toInt(strCallCode);
                    iCNState = Convert::toInt(strCNState);
                    add(strCallName, iCallCode, iCNState);
                }
            }
        }
    }

    file.close();

    // save index file name
    setIndexFile(strIndexFile);

    return true;
}

bool CallTable::add(CallRow& row)
{
    try {
        push_back(row);
    }
    catch (...) {
        return false;
    }
    return true;
}

bool CallTable::add(string strCallName, int iCallCode, int iCNState)
{
    try {
        CallRow row;
        row.setRowIndex(nextRowIndex());
        row.setCallName(strCallName);
        row.setCallCode(iCallCode);
        row.setCNState(iCNState);
        add(row);
    }
    catch (...) {
        return false;
    }
    return true;
}

void CallTable::sortByRowIndex()
{
    std::sort(begin(), end(), CallRow::compareRowIndex);
}

void CallTable::sortByCallName()
{ 
    std::sort(begin(), end(), CallRow::compareCallName);
}

void CallTable::sortByCallNameCNState()
{
    std::sort(begin(), end(), CallRow::compareCallNameCNState);
}

CallTable::iterator CallTable::find(std::string strCallName)
{
    CallTable::iterator iterFind;
    CallRow search;
    search.setCallName(strCallName);
    iterFind = lower_bound(begin(), end(), search, CallRow::compareCallName);
    if (iterFind != end() && iterFind->getCallName() == strCallName) {
        return iterFind;
    }
    else {
        return end();
    }
}

CallTable::iterator CallTable::find(std::string strCallName, int iCNState)
{
    CallTable::iterator iterFind;
    CallRow search;
    search.setCallName(strCallName);
    search.setCNState(iCNState);
    iterFind = lower_bound(begin(), end(), search, CallRow::compareCallName);
    if (iterFind != end() && iterFind->getCallName() == strCallName) {
        return iterFind;
    }
    else {
        return end();
    }
}

CallRow::CallRow()
{
    reset();
}

CallRow::~CallRow()
{
}

void CallRow::reset()
{
    m_iRowIndex = 0;
    m_strCallName = "";
    m_iCallCode = 0;
    m_iCNState = 0;
}

bool CallRow::compareRowIndex(const CallRow& elem1, const CallRow& elem2)
{
    return (elem1.m_iRowIndex < elem2.m_iRowIndex);
}

bool CallRow::compareCallName(const CallRow& elem1, const CallRow& elem2)
{
    return (elem1.m_strCallName < elem2.m_strCallName);
}

bool CallRow::compareCallNameCNState(const CallRow& elem1, const CallRow& elem2)
{
    if (elem1.m_strCallName != elem2.m_strCallName) { 
        return (elem1.m_strCallName < elem2.m_strCallName); 
    }
    return (elem1.m_iCNState < elem2.m_iCNState);
}

string CallRow::getRowStr()
{
    string strRow = "";
    strRow += "RowIndex: " + Convert::toString(m_iRowIndex);
    strRow += ", CallName: " + m_strCallName;
    strRow += ", CallCode: " + Convert::toString(m_iCallCode);
    strRow += ", CNState: " + Convert::toString(m_iCNState);    
    return strRow;
}
