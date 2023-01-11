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
// ~/apt2-src/util-package/IndexWriter.h ---
//
// Revision History:
// Created by David Le on 12/02/2015
//

#ifndef __INDEXWRITER_H__
#define __INDEXWRITER_H__

#include "IndexReader.h"

#include <time.h>
#include <iostream>
#include <string>
#include <algorithm>

enum OpenMode {
    OPENMODE_NONE = 0,
    OPENMODE_CREATE = 1,
    OPENMODE_UPDATE = 2,
};

class IndexWriter
{
public:
    IndexWriter();
    IndexWriter(std::string strIndexFile);
    IndexWriter(IndexReader* pindexReader);
    ~IndexWriter();

public:
    std::string getErrMsg() { return m_strErrMsg; }
    void setErrMsg(std::string str) { m_strErrMsg = str; }
    void clearErrMsg() { m_strErrMsg = ""; }

    std::string getIndexFile() { return m_strIndexFile; }
    void setIndexFile(std::string str) { m_strIndexFile = str; }

    FILE* getFileHandle() { return m_fileHandle; }
    void setFileHandle(FILE*& fp) { m_fileHandle = fp; }

    __int64 getFileSize() { return m_iFileSize; }
    void setFileSize(__int64 i) { m_iFileSize = i; }

    __int64 getFilePosition() { return m_iFilePosition; }
    void setFilePosition(__int64 i) { m_iFilePosition = i; }
    void incrementFilePosition(__int64 numBytes = 1) { m_iFilePosition += numBytes; }

    int getSampleMultiplier() { return m_iSampleMultiplier; }
    void setSampleMultiplier(int i) { m_iSampleMultiplier = i; }

    OpenMode getOpenMode() { return m_eOpenMode; }
    void setOpenMode(OpenMode e) { m_eOpenMode = e; }

    bool isOpen() { return m_bOpen; }
    void setOpen(bool b) { m_bOpen = b; }

    bool isLocalReader() { return m_bLocalReader; }
    void setLocalReader(bool b) { m_bLocalReader = b; }

    bool isOpenModeCreate() { return (m_eOpenMode == OPENMODE_CREATE); }
    bool isOpenModeUpdate() { return (m_eOpenMode == OPENMODE_UPDATE); }

    int getNumIndexRows() { return getIndexTable().getCount(); }

    IndexReader* getIndexReader() { return m_pindexReader; }
    IndexHeader& getIndexHeader() { return m_pindexReader->getIndexHeader(); }
    IndexTable& getIndexTable() { return m_pindexReader->getIndexTable(); }
    CallTable& getCallTable() { return getIndexHeader().getCallTable(); }

public:
    void reset();
    bool open();
    bool open(std::string strIndexFile);
    bool close();

    // blob index section

    bool createFile(std::string strIndexFile);

    bool write(const char* szData, size_t numBytes);
    bool writeStr(std::string line);
    bool writeLine(std::string line);

    bool writeFile();
    bool writeHeader();
    bool writeHeader_v1();
    bool writeBody();

    // standard operations

    double getFormatVersion() { return getIndexHeader().getFormatVersion(); }
    void setFormatVersion(double d) { getIndexHeader().setFormatVersion(d); }

    int getNumProbeSets() { return getIndexHeader().getNumProbeSets(); }
    void setNumProbeSets(int i) { getIndexHeader().setNumProbeSets(i); }
    void incrementNumProbeSets() { getIndexHeader().incrementNumProbeSets(); }

    int getNumSamples() { return getIndexHeader().getNumSamples(); }
    void setNumSamples(int i) { getIndexHeader().setNumSamples(i); }
    void incrementNumSamples() { getIndexHeader().incrementNumSamples(); }

    int getMaxAlleles() { return getIndexHeader().getMaxAlleles(); }
    void setMaxAlleles(int i) { getIndexHeader().setMaxAlleles(i); }
    
    int getMaxCNStates() { return getIndexHeader().getMaxCNStates(); }
    void setMaxCNStates(int i) { getIndexHeader().setMaxCNStates(i); }

    int getNumFields() { return getIndexHeader().getNumFields(); }
    void setNumFields(int i) { return getIndexHeader().setNumFields(i); }

    void clearFields() { getIndexHeader().getFields().clear(); }
    bool initFields(int maxAlleles);
    bool addField(std::string field, int size);
    void setFields(std::vector<std::pair<std::string, int> >& vFields) { getIndexHeader().getFields() = vFields; }

    void clearSampleNames() { getIndexHeader().getSampleNames().clear(); }
    bool initSampleNames(int numSamples);
    bool addSampleName(std::string strSampleName);

    void clearCallCodes() { getIndexHeader().getCallTable().clear(); }
    bool initCallCodes(int numCallCodes);
    bool addCallCode(std::string strCallName, int iCallCode, int iCNState);
    bool generateCallCodes();
    bool generateCallCodes_v1();

    void clearProbeSets() { getIndexTable().clear(); }
    bool initProbeSets(int numProbeSets);
    bool addProbeSet(std::string strProbeSetId, int numAlleles, __int64 fileOffset, int rowIndex = -1); // -1 for auto-increment

private:
    std::string m_strErrMsg;
    std::string m_strIndexFile;

    FILE* m_fileHandle;
    __int64 m_iFileSize;
    __int64 m_iFilePosition;
    int m_iSampleMultiplier;
    OpenMode m_eOpenMode;
    bool m_bOpen;
    bool m_bLocalReader;

    IndexReader* m_pindexReader;
};

#endif
