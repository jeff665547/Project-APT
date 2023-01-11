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
// ~/apt2-src/util-package/BlockReader.cpp ---
//
// Revision History:
// Created by David Le on 2/14/2016
//

#include "BlockReader.h"
#include "CodeMap.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include <iostream>
#include <sstream>
#include <string.h>
#include <memory>

using namespace std;

BlockReader::BlockReader()
{
    reset();
}

BlockReader::BlockReader(int iBlockSize)
{
    reset();
    m_iWindowSize = iBlockSize;
}

BlockReader::~BlockReader()
{
    close();
}

void BlockReader::reset()
{
    m_iWindowSize = 50;   // # of samples/columns to scan horizontally per window
    m_iWindowOffset = -1; // not loaded first sample in window
    m_iSampleIndex = -1;  // not loaded, first sample in blob

    m_vProbeSets.clear();
}

bool BlockReader::isReadNextWindow(int iSampleIndex)
{
    if (m_iWindowOffset == -1) {
        return true; // first window not loaded
    }
    if (iSampleIndex >= (m_iWindowOffset + m_iWindowSize)) {
        return true; // exceeds window limit, load next window
    }
    return false;
}

bool BlockReader::readNextWindow()
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to read data for next window, reader not opened.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to read data for next window, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (!seekToFirstProbeSet()) {
        setErrMsg("Failed to read data for next windw, unable to seek to first probeset.\tBlobFile: " + getBlobFile() + ", ErrMsg: " + getErrMsg());
        return false;
    }

    int iNextWindowOffset = (m_iWindowOffset < 0 ? 0 : m_iWindowOffset + m_iWindowSize);
    IndexTable& indexTable = m_indexReader.getIndexTable();
    m_vProbeSets.clear();
    m_vProbeSets.reserve(indexTable.getCount());
    for (int i = 0; i < indexTable.getCount(); i++) {
        IndexRow& row = indexTable.at(i);
        SampleList vSamples;
        if (!readNextWindowBlock(vSamples, row.getNumAlleles(), iNextWindowOffset, m_iWindowSize)) {
            setErrMsg("Failed to read data for next windw, unable to read next sample block.\tBlobFile: " + getBlobFile() + ", ErrMsg: " + getErrMsg());
            return false;
        }
        m_vProbeSets.push_back(make_pair(row, vSamples));
    }

    // move offset to next window of samples
    m_iWindowOffset = iNextWindowOffset;

    return true;
}

bool BlockReader::nextSample()
{
    if (!isOpen()) {
        return false;
    }

    m_iSampleIndex++;
    if (0 <= m_iSampleIndex && m_iSampleIndex < getNumSamples()) {
        if (isReadNextWindow(m_iSampleIndex)) {
            return readNextWindow();
        }
        return true;
    }

    return false;
}

string BlockReader::getSampleName()
{ 
    if (0 <= m_iSampleIndex && m_iSampleIndex < getNumSamples()) {
        return getSampleNames()[m_iSampleIndex]; 
    }
    return "";
}

string BlockReader::getSampleCall(SampleList& vSamples)
{
    bool bSampleInBlobRange = (0 <= m_iSampleIndex && m_iSampleIndex < getNumSamples());
    bool bSampleInWindowRange = (m_iWindowOffset <= m_iSampleIndex && m_iSampleIndex <= (m_iWindowOffset + m_iWindowSize));
    if (bSampleInBlobRange && bSampleInWindowRange) {
        int iWindowIndex = (m_iSampleIndex % m_iWindowSize);
        if (0 <= iWindowIndex && iWindowIndex < vSamples.getCount()) {
            uint8_t byCallCode = vSamples[iWindowIndex].getCall();
            //
            // todo: should use CallCode table in IndexHeader section instead
            //
            switch (byCallCode) {
                case 6: return "A/A";
                case 7: return "B/B";
                case 8: return "A/B";
                case 11: return "NoCall";
                case 12: return "OTV";
                default: return "NoCall";
            }
        }
        else {
            setErrMsg("Failed to extract call from sample for probeset, window index is out of range.\tBlobFile: " + getBlobFile() + ", WindowIndex: " + Convert::toString(iWindowIndex) + ", SamplesInWindow: " + Convert::toString(vSamples.getCount()));
            return "NoCall";
        }
    }
    else {
        setErrMsg("Failed to extract call from sample for probeset, sample index is out of range.\tBlobFile: " + getBlobFile() + ", SampleIndex: " + Convert::toString(m_iSampleIndex) + ", NumSamples: " + Convert::toString(getNumSamples()) + ", WindowOffset: " + Convert::toString(m_iWindowOffset) + ", WindowSize: " + Convert::toString(m_iWindowSize));
        return "NoCall";
    }
}
