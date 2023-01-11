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
// ~/apt2-src/util-package/SampleTypesEx.h ---
//
// Revision History:
// Created by David Le on 11/20/2015
//

#ifndef __SAMPLETYPESEX_H__
#define __SAMPLETYPESEX_H__

#include "portability/affy-base-types.h"

#include <cstring>
#include <string>
#include <string.h>
#include <vector>

// keep actual size of packed structure
#pragma pack(push, 1)
class SampleEx
{
public:
    SampleEx(int iNumAlleles = 2);
	~SampleEx();

public:
    uint8_t getCall() { return m_byCall; }
    void setCall(uint8_t by) { m_byCall = by; }
    int getCallInt() { return (int8_t)m_byCall; }

	float getConfidence() { return m_fConfidence; }
	void setConfidence(float f) { m_fConfidence = f; }

	float getLogRatio() { return m_fLogRatio; }
	void setLogRatio(float f) { m_fLogRatio = f; }
	
	float getStrength() { return m_fStrength; }
	void setStrength(float f) { m_fStrength = f; }

    int getNumAlleles() { return (int)m_vAlleleSignals.size(); }
    std::vector<float>& getAlleleSignals() { return m_vAlleleSignals; }
    float getAlleleSignal(int index) { return m_vAlleleSignals[index]; }
    void setAlleleSignal(int index, float signal) { m_vAlleleSignals[index] = signal; }

public:
	void reset();
    void init(int iNumAlleles);
    size_t format(char* szData, int iSampleMultiplier = 1);
    size_t formatNext(char*& szData, int iSampleMultiplier = 1);
    bool parse(char* szData, int numAlleles);
    bool parseNext(char*& szData, int numAlleles);
    float computeLogRatio(); // contrast
    float computeStrength(); // size
    std::string getSampleStr();

private:
    uint8_t m_byCall;
    float m_fConfidence;
    float m_fLogRatio;
    float m_fStrength;
    std::vector<float> m_vAlleleSignals;
};
#pragma pack(pop)

class SampleList : public std::vector<SampleEx>
{
public:
    SampleList();
    SampleList(int iNumSamples, int iNumAlleles = 2);
    ~SampleList();

public:
    int getCount() { return (int)size(); }
    bool formatBlock(char* szData, int iSampleMultiplier = 1);
    bool parseBlock(char* szData, int numSamples, int numAlleles);
    bool parseWindow(char* szData, int numSamples, int numAlleles, int bytesPerSample, int windowOffset, int windowSize);
    std::string getSamplesStr();
};

class CallEx
{
public:
    CallEx();
    CallEx(uint8_t byCall);
    ~CallEx();

public:
    uint8_t getCall() { return m_byCall; }
    void setCall(uint8_t by) { m_byCall = by; }

public:
    void reset();

private:
    uint8_t m_byCall;
};

class CallList : public std::vector<CallEx>
{
public:
    CallList();
    ~CallList();

public:
    size_t getCount() { return size(); }
};

class AlleleSignalEx
{
public:
    AlleleSignalEx(int iNumAlleles = 2);
    AlleleSignalEx(float fAAlleleSignal, float fBAlleleSignal);
    ~AlleleSignalEx();

public:
    float getAlleleSignal(int index) { return (float)m_arAlleleSignals[index * sizeof(float)]; }
    void setAlleleSignal(int index, float f) { memcpy(&m_arAlleleSignals[index * sizeof(float)], &f, sizeof(float)); }

public:
    void reset();
    void init(int iNumAlleles);

private:
    uint8_t* m_arAlleleSignals;
    int m_iNumAlleles;
};

class AlleleSignalList : public std::vector<AlleleSignalEx>
{
public:
    AlleleSignalList(int iNumAlleles = 2);
    AlleleSignalList(int iNumSamples, int iNumAlleles = 2);
    ~AlleleSignalList();

public:
    int getNumAlleles() { return m_iNumAlleles; }
    void setNumAlleles(int i) { m_iNumAlleles = i; }

public:
    size_t getCount() { return size(); }

private:
    int m_iNumAlleles;
};

#endif
