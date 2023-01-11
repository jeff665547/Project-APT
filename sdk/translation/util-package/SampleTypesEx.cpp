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
// ~/apt2-src/util-package/SampleTypesEx.cpp ---
//
// Revision History:
// Created by David Le on 11/20/2015
//

#include "SampleTypesEx.h"
#include "util/Convert.h"
#include "util/Util.h"

#include <cstring>
#include <string.h>

using namespace std;

SampleEx::SampleEx(int iNumAlleles)
{
    reset();
    init(iNumAlleles);
}

SampleEx::~SampleEx()
{
}

void SampleEx::reset()
{
    m_byCall = (uint8_t)0;
    m_fConfidence = 0.0f;
    m_fLogRatio = 0.0f;
    m_fStrength = 0.0f;
    m_vAlleleSignals.clear();
}

void SampleEx::init(int iNumAlleles)
{
    m_vAlleleSignals.assign(iNumAlleles, 0.0f);
}

size_t SampleEx::format(char* szData, int iSampleMultiplier)
{
    char* ptr = szData;
    return formatNext(ptr, iSampleMultiplier);
}

size_t SampleEx::formatNext(char*& szData, int iSampleMultiplier)
{
    char* ptr = szData;
    try {
        int sizeByte = sizeof(uint8_t);
        int sizeFloat = sizeof(float);
        for (int i = 0; i < iSampleMultiplier; i++) {
            std::memcpy(szData, &m_byCall, sizeByte); szData += sizeByte;
            std::memcpy(szData, &m_fConfidence, sizeFloat); szData += sizeFloat;
            std::memcpy(szData, &m_fLogRatio, sizeFloat); szData += sizeFloat;
            std::memcpy(szData, &m_fStrength, sizeFloat); szData += sizeFloat;
            for (int j = 0; j < getNumAlleles(); j++) {
                std::memcpy(szData, &m_vAlleleSignals[j], sizeFloat);
                szData += sizeFloat;
            }
        }
    }
    catch (...) {
        return -1;
    }
    return (szData - ptr);
}

bool SampleEx::parse(char* szData, int numAlleles)
{
    char* ptr = szData;
    return parseNext(ptr, numAlleles);
}

bool SampleEx::parseNext(char*& szData, int numAlleles)
{
    try {
        init(numAlleles);
        int sizeByte = sizeof(uint8_t);
        int sizeFloat = sizeof(float);
        std::memcpy(&m_byCall, szData, sizeByte); szData += sizeByte;
        std::memcpy(&m_fConfidence, szData, sizeFloat); szData += sizeFloat;
        std::memcpy(&m_fLogRatio, szData, sizeFloat); szData += sizeFloat;
        std::memcpy(&m_fStrength, szData, sizeFloat); szData += sizeFloat;
        for (int i = 0; i < getNumAlleles(); i++) {
            std::memcpy(&m_vAlleleSignals[i], szData, sizeFloat);
            szData += sizeFloat;
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

float SampleEx::computeLogRatio() // contrast
{
    float fAAlleleSignal = getAlleleSignal(0);
    float fBAlleleSignal = getAlleleSignal(1);
    m_fLogRatio = (log2(fBAlleleSignal) - log2(fAAlleleSignal)) * -1.0f;
	return m_fLogRatio;
}

float SampleEx::computeStrength() // size
{
    float fAAlleleSignal = getAlleleSignal(0);
    float fBAlleleSignal = getAlleleSignal(1);
    m_fStrength = (log2(fAAlleleSignal) + log2(fBAlleleSignal)) / 2.0f;
	return m_fStrength;
}

string SampleEx::getSampleStr()
{
    string str = "";
    str += "Call: " + Convert::toString(m_byCall);
    str += ", Confidence: " + Convert::toString(m_fConfidence);
    str += ", LogRatio: " + Convert::toString(m_fLogRatio);
    str += ", Strength: " + Convert::toString(m_fStrength);
    for (int i = 0; i < getNumAlleles(); i++) {
        str += ", ";
        str += ('A' + i); 
        str += "AlleleSignal: " + Convert::toString(getAlleleSignal(i));
    }
    return str;
}

SampleList::SampleList()
{
}

SampleList::SampleList(int iNumSamples, int iNumAlleles) 
{
    assign(iNumSamples, SampleEx(iNumAlleles));
}

SampleList::~SampleList()
{
}

bool SampleList::formatBlock(char* szData, int iSampleMultiplier)
{
    try {
        char* ptr = szData;
        for (int i = 0; i < getCount(); i++) {
            SampleEx& sample = at(i);
            if (!sample.formatNext(ptr, iSampleMultiplier)) {
                return false;
            }
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

bool SampleList::parseBlock(char* szData, int numSamples, int numAlleles)
{
    try {
        clear();
        reserve(numSamples);
        char* ptr = szData;
        for (int i = 0; i < numSamples; i++) {
            SampleEx sample;
            if (!sample.parseNext(ptr, numAlleles)) {
                return false;
            }
            push_back(sample);
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

bool SampleList::parseWindow(char* szData, int numSamples, int numAlleles, int bytesPerSample, int windowOffset, int windowSize)
{
    try {
        clear();
        reserve(windowSize);
        char* ptr = szData;
        ptr += (windowOffset * bytesPerSample);
        for (int i = windowOffset; i < (windowOffset + windowSize) && i < numSamples; i++) {
            SampleEx sample;
            if (!sample.parseNext(ptr, numAlleles)) {
                return false;
            }
            push_back(sample);
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

string SampleList::getSamplesStr()
{
    string str = "[";
    for (int i = 0; i < getCount(); i++) {
        string strSample = at(i).getSampleStr();
        str += "{" + strSample + "},";
    }
    if (str.at(str.length() - 1) == ',') {
        str = str.substr(0, str.length() - 1);
    }
    str += "]";
    return str;
}

CallEx::CallEx()
{
    reset();
}

CallEx::CallEx(uint8_t byCall)
{
    reset();
    m_byCall = byCall;
}

CallEx::~CallEx()
{
}

void CallEx::reset()
{
    m_byCall = (uint8_t)0;
}

CallList::CallList()
{
}

CallList::~CallList()
{
}

AlleleSignalEx::AlleleSignalEx(int iNumAlleles)
{
    m_arAlleleSignals = NULL;
    init(iNumAlleles);
}

AlleleSignalEx::AlleleSignalEx(float fAAlleleSignal, float fBAlleleSignal)
{
    m_arAlleleSignals = NULL;
    init(2);
    setAlleleSignal(0, fAAlleleSignal);
    setAlleleSignal(1, fBAlleleSignal);
}

AlleleSignalEx::~AlleleSignalEx()
{
    reset();
}

void AlleleSignalEx::reset()
{
    if (m_arAlleleSignals != NULL) {
        delete[] m_arAlleleSignals;
        m_arAlleleSignals = NULL;
    }
    m_iNumAlleles = 0;
}

void AlleleSignalEx::init(int iNumAlleles)
{
    reset();
    m_arAlleleSignals = new uint8_t[iNumAlleles * sizeof(float)];
    m_iNumAlleles = iNumAlleles;
}

AlleleSignalList::AlleleSignalList(int iNumAlleles)
{
    setNumAlleles(iNumAlleles);
}

AlleleSignalList::AlleleSignalList(int iNumSamples, int iNumAlleles)
{
    assign(iNumSamples, AlleleSignalEx(iNumAlleles));
    setNumAlleles(iNumAlleles);
}

AlleleSignalList::~AlleleSignalList()
{
}
