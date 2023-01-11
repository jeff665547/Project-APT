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
// ~/apt2-src/util-package/BlockReader.h ---
//
// Revision History:
// Created by David Le on 2/14/2016
//

#ifndef __BLOCKREADER_H__
#define __BLOCKREADER_H__

#include "BlobReader.h"
#include "SampleTypesEx.h"
#include <map>
#include <set>
#include <vector>

class BlockReader : public BlobReader
{
public:
    BlockReader();
    BlockReader(int iBlockSize);
    ~BlockReader();

public:
    int getWindowSize() { return m_iWindowSize; }
    void setWindowSize(int i) { m_iWindowSize = i; }

    int getWindowOffset() { return m_iWindowOffset; }
    void setWindowOffset(int i) { m_iWindowOffset = i; }

    int getSampleIndex() { return m_iSampleIndex; }
    void setSampleIndex(int i) { m_iSampleIndex = i; }

    std::vector<std::pair<IndexRow, SampleList> >& getProbeSets() { return m_vProbeSets; }

public:
    void reset();
    bool isReadNextWindow(int iSampleIndex);
    bool readNextWindow();
    bool nextSample();
    std::string getSampleName();
    std::string getSampleCall(SampleList& vSamples);

private:
    int m_iWindowSize;   // # of samples to scan horizontally per window
    int m_iWindowOffset; // index of first sample in window
    int m_iSampleIndex;

    std::vector<std::pair<IndexRow, SampleList> > m_vProbeSets;
};

#endif
