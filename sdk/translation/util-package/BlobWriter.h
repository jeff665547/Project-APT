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
// ~/apt2-src/util-package/BlobWriter.h ---
//
// Revision History:
// Created by David Le on 12/02/2015
//

#ifndef __BLOBWRITER_H__
#define __BLOBWRITER_H__

#include "BlobReader.h"
#include "IndexWriter.h"
#include "PackageConsts.h"
#include <time.h>
#include <iostream>
#include <string>
#include <algorithm>

class BlobWriter
{
public:
    static size_t getFileSize(std::string strFile) { return BlobReader::getFileSize(strFile); }

public:
    BlobWriter();
    BlobWriter(int maxAlleles);
    ~BlobWriter();

public:
    std::string getErrMsg() { return m_strErrMsg; }
    void setErrMsg(std::string str) { m_strErrMsg = str; }
    void clearErrMsg() { m_strErrMsg = ""; }

    std::string getBlobFile() { return m_strBlobFile; }
    void setBlobFile(std::string str) { m_strBlobFile = str; }

    std::string getIndexFile() { return m_strIndexFile; }
    void setIndexFile(std::string str) { m_strIndexFile = str; }

    std::string getBatchDataDir() { return m_strBatchDataDir; }
    void setBatchDataDir(std::string str) { m_strBatchDataDir = str; }

    FILE* getFileHandle() { return m_fileHandle; }
    void setFileHandle(FILE*& fp) { m_fileHandle = fp; }

    __int64 getFileSize() { return m_iFileSize; }
    void setFileSize(__int64 i) { m_iFileSize = i; }

    __int64 getFilePosition() { return m_iFilePosition; }
    void setFilePosition(__int64 i) { m_iFilePosition = i; }
    void incrementFilePosition(__int64 numBytes = 1) { m_iFilePosition += numBytes; }

    OpenMode getOpenMode() { return m_eOpenMode; }
    void setOpenMode(OpenMode e) { m_eOpenMode = e; }

    bool isOpen() { return m_bOpen; }
    void setOpen(bool b) { m_bOpen = b; }

    __int64 getWrittenSampleBlocks() { return m_iWrittenSampleBlocks; }
    void setWrittenSampleBlocks(__int64 i) { m_iWrittenSampleBlocks = i; }
    void incrementWrittenSampleBlocks(__int64 numSamples = 1) { m_iWrittenSampleBlocks += numSamples; }

    __int64 getWrittenNumAlleles() { return m_iWrittenNumAlleles; }
    void setWrittenNumAlleles(__int64 i) { m_iWrittenNumAlleles = i; }
    void incrementWrittenNumAlleles(__int64 numAlleles = 1) { m_iWrittenNumAlleles += numAlleles; }
    
    char* getBuffer() { return m_szBuffer; }
    void setBuffer(char* data, int len) { std::memcpy(m_szBuffer, data, len); }

    int getBufferLength() { return m_iBufferLength; }
    void setBufferLength(int i) { m_iBufferLength = i; }

    int getBufferPosition() { return m_iBufferPosition; }
    void setBufferPosition(int i) { m_iBufferPosition = i; }
    void incrementBufferPosition(int numBytes = 1) { m_iBufferPosition += numBytes; }
    
    bool isEndOfFile() { return m_bEndOfFile; }
    void setEndOfFile(bool b) { m_bEndOfFile = b; }

    int getSampleMultiplier() { return m_indexWriter.getSampleMultiplier(); }
    void setSampleMultiplier(int i) { m_indexWriter.setSampleMultiplier(i); }

    bool isOpenModeCreate() { return (m_eOpenMode == OPENMODE_CREATE); }
    bool isOpenModeUpdate() { return (m_eOpenMode == OPENMODE_UPDATE); }

    BlobReader& getBlobReader() { return m_blobReader; }
    IndexWriter& getIndexWriter() { return m_indexWriter; }
    IndexHeader& getIndexHeader() { return m_indexWriter.getIndexHeader(); }
    CallTable& getCallTable() { return m_indexWriter.getCallTable(); }
    IndexTable& getIndexTable() { return m_indexWriter.getIndexTable(); }

public:
    void reset();
    void resetBuffer();
    bool initFolder(std::string strPackageFolder);
    bool init(std::string strBlobFile, std::string strIndexFile = "");
    bool open();
    bool open(std::string strBlobFile, std::string strIndexFile);
    bool close();

    bool write(const char* szData, size_t numBytes);
    bool writeByte(uint8_t byValue) { return write((char*)&byValue, sizeof(uint8_t)); }
    bool writeFloat(float fValue) { return write((char*)&fValue, sizeof(float)); }
    bool writeStr(std::string line);
    bool writeLine(std::string line);
    bool writeBuffer(const char* szData, int numBytes);
    bool writeBufferStr(std::string str);
    bool writeBufferLine(std::string line);
    bool flushBuffer();

    // blob data/index section

    size_t getCurrentFilePosition();
    bool seek(size_t fileOffset);
    bool seekToFirstProbeSet();
    bool seekToProbeSet(std::string strProbeSetId);
    bool seekToSample(std::string strProbeSetId, std::string strSampleName);
    bool skipToNextSample(int numAlleles, int numSamples = 1);
    bool skipToNextSampleBlock(int numAlleles);
    bool writeSample(SampleEx& sample);
    bool writeSampleBlock(SampleList& vSamples);
    bool writeNextSample(SampleEx& sample);
    bool writeNextSampleBlock(SampleList& vSamples);
    
    // standard operations

    double getFormatVersion() { return getIndexHeader().getFormatVersion(); }
    void setFormatVersion(double d) { getIndexHeader().setFormatVersion(d); }

    std::string getTransformMethod() { return getIndexHeader().getTransformMethod(); }
    void setTransformMethod(std::string str) { getIndexHeader().setTransformMethod(str); }

    int getNumProbeSets() { return m_indexWriter.getNumProbeSets(); }
    void setNumProbeSets(int i) { m_indexWriter.setNumProbeSets(i); }

    int getNumSamples() { return m_indexWriter.getNumSamples(); }
    void setNumSamples(int i) { m_indexWriter.setNumSamples(i); }
    int getNumRawSamples() { return (getNumSamples() / getSampleMultiplier()); }

    int getMaxAlleles() { return m_indexWriter.getMaxAlleles(); }
    void setMaxAlleles(int i) { m_indexWriter.setMaxAlleles(i); }

    int getMaxCNStates() { return m_indexWriter.getMaxCNStates(); }
    void setMaxCNStates(int i) { m_indexWriter.setMaxCNStates(i); }

    int getNumFields() { return getIndexHeader().getNumFields(); }
    int getNumCallRows() { return getCallTable().getCount(); }
    int getNumIndexRows() { return getIndexTable().getCount(); }

    bool initSampleNames(int numSamples) { return m_indexWriter.initSampleNames(numSamples); }
    bool addSampleName(std::string strSampleName) { return m_indexWriter.addSampleName(strSampleName); }

    bool initCallCodes(int numCallCodes) { return m_indexWriter.initCallCodes(numCallCodes); }
    bool addCallCode(std::string strCallName, int iCallCode, int iCNState) { return m_indexWriter.addCallCode(strCallName, iCallCode, iCNState); }
    bool generateCallCodes() { return m_indexWriter.generateCallCodes(); }

    bool initProbeSets(int numProbeSets) { return m_indexWriter.initProbeSets(numProbeSets); }
    bool addProbeSet(std::string strProbeSetId, int numAlleles) { return m_indexWriter.addProbeSet(strProbeSetId, numAlleles, m_iFilePosition); }

    bool initFields(int maxAlleles) { return m_indexWriter.initFields(maxAlleles); }

    __int64 getExpectedSampleBlocks() { return (__int64)getNumProbeSets() * (__int64)getNumSamples(); }
    __int64 getExpectedFileSize() { return (getExpectedSampleBlocks() * 13) + (getWrittenNumAlleles() * 4); }
    __int64 getTotalNumAlleles() { return (getWrittenNumAlleles() / getNumSamples()); }

private:
    std::string m_strErrMsg;
    std::string m_strBlobFile;
    std::string m_strIndexFile;
    std::string m_strBatchDataDir;

    FILE* m_fileHandle;
    __int64 m_iFileSize;
    __int64 m_iFilePosition;
    OpenMode m_eOpenMode;
    bool m_bOpen;
    __int64 m_iWrittenSampleBlocks;
    __int64 m_iWrittenNumAlleles; // # of alleles written for all samples in all probesets

    char m_szBuffer[MAX_BUFFER];
    int m_iBufferLength;
    int m_iBufferPosition;
    bool m_bEndOfFile;

    BlobReader m_blobReader;
    IndexWriter m_indexWriter;
};

#endif
