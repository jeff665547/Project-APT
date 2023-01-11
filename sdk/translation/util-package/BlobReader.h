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
// ~/apt2-src/util-package/BlobReader.h ---
//
// Revision History:
// Created by David Le on 12/01/2015
//

#ifndef __BLOBREADER_H__
#define __BLOBREADER_H__

#include "IndexReader.h"
#include "PackageConsts.h"
#include "SampleTypesEx.h"
#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <set>
#include <vector>

class BlobReader
{
public :
    static char* readFile(std::string strBlobFile, size_t& numBytes);
    static size_t getFileSize(std::string strFile);
    static bool is64BitFile(std::string strFile);
    static int bytesForField(std::string strFieldType);

public:
    BlobReader();
    ~BlobReader();

public:
    std::string getErrMsg() { return m_strErrMsg; }
    void setErrMsg(std::string str) { m_strErrMsg = str; }
    void clearErrMsg() { m_strErrMsg = ""; }

    std::string getBlobFile() { return m_strBlobFile; }
    void setBlobFile(std::string str) { m_strBlobFile = str; }

    std::string getIndexFile() { return m_strIndexFile; }
    void setIndexFile(std::string str) { m_strIndexFile = str; }

    FILE* getFileHandle() { return m_fileHandle; }
    void setFileHandle(FILE*& fp) { m_fileHandle = fp; }

    __int64 getFileSize() { return m_iFileSize; }
    void setFileSize(__int64 i) { m_iFileSize = i; }

    __int64 getFilePosition() { return m_iFilePosition; }
    void setFilePosition(__int64 i) { m_iFilePosition = i; }
    void incrementFilePosition(__int64 numBytes = 1) { m_iFilePosition += numBytes; }

    bool isEndOfFile() { return m_bEndOfFile; }
    void setEndOfFile(bool b) { m_bEndOfFile = b; }

    bool isOpen() { return m_bOpen; }
    void setOpen(bool b) { m_bOpen = b; }

    char* getBuffer() { return m_szBuffer; }
    void setBuffer(char* data, int len) { std::memcpy(m_szBuffer, data, len); }

    int getBufferLength() { return m_iBufferLength; }
    void setBufferLength(int i) { m_iBufferLength = i; }

    int getBufferPosition() { return m_iBufferPosition; }
    void setBufferPosition(int i) { m_iBufferPosition = i; }
    void incrementBufferPosition(int numBytes = 1) { m_iBufferPosition += numBytes; }

    int getNumCallRows() { return getIndexHeader().getNumCallRows(); }

    IndexReader& getIndexReader() { return m_indexReader; }
    IndexHeader& getIndexHeader() { return m_indexReader.getIndexHeader(); }
    IndexTable& getIndexTable() { return m_indexReader.getIndexTable(); }
    CallTable& getCallTable() { return getIndexHeader().getCallTable(); }

public:
    void reset();
    void resetBuffer();
    bool initFolder(std::string strPackageFolder);
    bool init(std::string strBlobFile, std::string strIndexFile = "");
    bool open();
    bool open(std::string strBlobFile, std::string strIndexFile = "");
    bool close();
    
    std::string readLine();
    std::string nextBufferLine();
    bool loadBufferData();

    // blob data/index section

    size_t getCurrentFilePosition();
    bool seek(size_t fileOffset);
    bool seekToFirstProbeSet();
    bool seekToProbeSet(std::string strProbeSetId);
    bool seekToSample(std::string strProbeSetId, std::string strSampleName);
    bool skipToNextSample(int numAlleles, int numSamples = 1);
    bool skipToNextSampleBlock(int numAlleles);
    bool readNextSample(SampleEx& sample, int numAlleles);
    bool readNextSampleBlock(SampleList& vSamples, int numAlleles);
    bool readNextWindowBlock(SampleList& vSamples, int numAlleles, int windowOffset, int windowSize);
    bool readIndexFile(std::string strIndexFile);

    // standard operations

    double getFormatVersion() { return getIndexHeader().getFormatVersion(); }
    void setFormatVersion(double d) { getIndexHeader().setFormatVersion(d); }
    std::string getFormatVersionStr() { char sz[MAX_LINE]; sprintf(sz, "%.1lf", getFormatVersion()); return sz; }

    std::string getTransformMethod() { return getIndexHeader().getTransformMethod(); }
    void setTransformMethod(std::string str) { getIndexHeader().setTransformMethod(str); }

    int getNumProbeSets() { return getIndexHeader().getNumProbeSets(); }
    void setNumProbeSets(int i) { getIndexHeader().setNumProbeSets(i); }

    int getNumSamples() { return getIndexHeader().getNumSamples(); }
    void setNumSamples(int i) { getIndexHeader().setNumSamples(i); }

    int getNumFields() { return getIndexHeader().getNumFields(); }
    void setNumFields(int i) { getIndexHeader().setNumFields(i); }

    int getMaxAlleles() { return getCallTable().getMaxAlleles(); }
    void setMaxAlleles(int i) { return getCallTable().setMaxAlleles(i); }

    int getMaxCNStates() { return getCallTable().getMaxCNStates(); }
    void setMaxCNStates(int i) { return getCallTable().setMaxCNStates(i); }

    int getNumIndexRows() { return getIndexTable().getCount(); }
    int getNumCallCodes() { return getCallTable().getCount(); }

    std::string getFieldName(int index) { return getFields().at(index).first; }
    std::string getCallCodeStr(int index) { return getCallTable().at(index).getCallCodeStr(); }
    std::string getCallName(int index) { return getCallTable().at(index).getCallName(); }
    int getCallCode(int index) { return getCallTable().at(index).getCallCode(); }
    int getCNState(int index) { return getCallTable().at(index).getCNState(); }

    int getBytesPerSample(int numAlleles) { return getIndexHeader().getBytesPerSample(numAlleles); }
    int getBytesPerSampleBlock(int numAlleles) { return getIndexHeader().getBytesPerSampleBlock(numAlleles); }
    int getBytesPerSampleBlock(SampleList& vSamples) { return getIndexHeader().getBytesPerSampleBlock(vSamples); }
    int findSampleIndex(std::string strSampleName) { return getIndexHeader().findSampleIndex(strSampleName); }

    std::vector<std::string>& getSampleNames() { return getIndexHeader().getSampleNames(); }
    std::vector<std::pair<std::string, int> >& getFields() { return getIndexHeader().getFields(); }
    std::map<std::string, int>& getMapSampleToIndexes() { return getIndexHeader().getMapSampleToIndexes(); }
    std::map<std::string, int>& getMapFieldToIndexes() { return getIndexHeader().getMapFieldToIndexes(); }

protected:
    std::string m_strErrMsg;
    std::string m_strBlobFile;
    std::string m_strIndexFile;

    FILE* m_fileHandle;
    __int64 m_iFileSize;
    __int64 m_iFilePosition;
    bool m_bOpen;

    char m_szBuffer[MAX_BUFFER];
    int m_iBufferLength;
    int m_iBufferPosition;
    bool m_bEndOfFile;

    IndexReader m_indexReader;
};

#endif
