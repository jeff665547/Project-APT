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
// ~/apt2-src/util-package/IndexWriter.cpp ---
//
// Revision History:
// Created by David Le on 12/04/2015
//

#include "IndexWriter.h"
#include "util/Convert.h"
#include "util/Fs.h"

#include <iostream>
#include <string>

using namespace std;

#ifdef __APPLE__
#define fseeko64 fseeko
#define ftello64 ftello
#endif

IndexWriter::IndexWriter()
{
    m_pindexReader = NULL;
    reset();
    m_pindexReader = new IndexReader();
    setLocalReader(true);
}

IndexWriter::IndexWriter(string strIndexFile)
{
    m_pindexReader = NULL;
    reset();
    m_pindexReader = new IndexReader();
    setLocalReader(true);
    setIndexFile(strIndexFile);
}

IndexWriter::IndexWriter(IndexReader* pindexReader)
{
    m_pindexReader = NULL;
    reset();
    m_pindexReader = pindexReader;
    setLocalReader(false);
    setIndexFile(pindexReader->getIndexFile());
}

IndexWriter::~IndexWriter()
{
    if (isLocalReader()) {
        delete m_pindexReader;
        m_pindexReader = NULL;
    }
    close();
}

void IndexWriter::reset()
{
    m_strErrMsg = "";
    m_strIndexFile = "";
    
    m_fileHandle = NULL;
    m_iFileSize = 0;
    m_iFilePosition = 0;
    m_iSampleMultiplier = 1;
    m_eOpenMode = OPENMODE_NONE;
    m_bOpen = false;
    m_bLocalReader = false;

    if (m_pindexReader != NULL) {
        m_pindexReader->reset();
    }
}

bool IndexWriter::open()
{
    return open(getIndexFile());
}

bool IndexWriter::open(string strIndexFile)
{
    clearErrMsg();

    try {
        setOpen(false);

        if (Fs::fileExists(strIndexFile)) {
            setOpenMode(OPENMODE_UPDATE);
            setFilePosition(0);
            m_fileHandle = fopen(strIndexFile.c_str(), "r+b"); // update index file
            if (!m_fileHandle) {
                setErrMsg("Failed to open index file for updating, file i/o error.\tIndexFile: " + strIndexFile);
                return false;
            }
            // parse index file
            if (isLocalReader() && m_pindexReader->readFile(strIndexFile)) {
                setErrMsg(m_pindexReader->getErrMsg());
                return false;
            }
        }
        else {
            setOpenMode(OPENMODE_CREATE);
            setFilePosition(0);
            m_fileHandle = fopen(strIndexFile.c_str(), "w+b"); // create index file
            if (!m_fileHandle) {
                setErrMsg("Failed to open index file for writing, file i/o error.\tIndexFile: " + strIndexFile);
                return false;
            }
        }
    }
    catch (...) {
        setErrMsg("Exception occurred while opening index writer for processing.\tIndexFile: " + strIndexFile);
        return false;
    }

    // index file is opened
    setOpen(true);
    setIndexFile(strIndexFile);

    return true;
}

bool IndexWriter::close()
{
    clearErrMsg();

    // close blob file handle
    if (m_fileHandle != NULL) {
        fclose(m_fileHandle);
        m_fileHandle = NULL;
    }

    reset();

    return true;
}

bool IndexWriter::createFile(string strIndexFile)
{
    clearErrMsg();

    try {
        if (!open(strIndexFile)) {
            return false;
        }

        if (!writeFile()) {
            return false;
        }

        if (!close()) {
            return false;
        }
    }
    catch (...) {
        return false;
    }

    return true;
}

bool IndexWriter::write(const char* szData, size_t numBytes)
{
    clearErrMsg();

    try {
        if (!isOpen()) {
            setErrMsg("Failed to write index to stream, writer not opened.\tIndexFile: " + getIndexFile());
            return false;
        }

        if (m_fileHandle == NULL) {
            setErrMsg("Failed to write index to stream, file handle is invalid.\tIndexFile: " + getIndexFile() + ", NumBytes: " + Convert::toString(numBytes));
            return false;
        }

        size_t bytesWritten = fwrite(szData, 1, numBytes, m_fileHandle);
        if (bytesWritten != numBytes) {
            setErrMsg("Failed to write data to index file, file i/o error.\tIndexFile: " + getIndexFile() + ", NumBytes: " + Convert::toString(numBytes));
            return false;
        }

        incrementFilePosition(numBytes);
    }
    catch (...) {
        setErrMsg("Exception occurred while writing index to stream.\tIndexFile: " + getIndexFile() + ", NumBytes: " + Convert::toString(numBytes));
        return false;
    }

    return true;
}

bool IndexWriter::writeStr(string str)
{
    return write(str.c_str(), str.size());
}

bool IndexWriter::writeLine(string line)
{
    line += "\n";
    return write(line.c_str(), line.size());
}

bool IndexWriter::writeFile()
{
    bool bSuccessful = true;

    clearErrMsg();

    try {
        if (!isOpen()) {
            setErrMsg("Failed to write content to index file, writer not opened.\tIndexFile: " + getIndexFile());
            return false;
        }

        if (m_fileHandle == NULL) {
            setErrMsg("Failed to write content to index file, file handle is invalid.\tIndexFile: " + getIndexFile());
            return false;
        }

        if (!writeHeader()) {
            bSuccessful = false;
        }

        if (!writeBody()) {
            bSuccessful = false;
        }
    }
    catch (...) {
        setErrMsg("Exception occurred while writing content to index file.\tIndexFile: " + getIndexFile());
        return false;
    }

    return bSuccessful;
}

// index header section
bool IndexWriter::writeHeader()
{
    if (getFormatVersion() == 1.0) {
        return writeHeader_v1();
    }

    bool bSuccessful = true;

    writeLine("#%%format-version=" + getIndexHeader().getFormatVersionStr());
    writeLine("#%%brlmmp-transform=" + getIndexHeader().getTransformMethod());
    writeLine("#");
    writeLine("#%%num-probesets=" + Convert::toString(getNumIndexRows()));
    writeLine("#");

    // write sample names
    vector<string>& vSampleNames = getIndexHeader().getSampleNames();
    int numSamples = (int)vSampleNames.size();
    writeLine("#%%num-samples=" + Convert::toString(numSamples));
    if (numSamples > 0) {
        int count = 1;
        for (int i = 0; i < numSamples; i++) {
            string strSampleName = vSampleNames[i];
            writeLine("#%%cel-name-" + Convert::toString(count) + "=" + strSampleName);
            count++;
        }
    }
    else {
        setErrMsg("Failed to write header section for index file, no samples detected.\tIndexFile: " + getIndexFile());
        bSuccessful = false;
    }
    writeLine("#");

    // write sample fields
    vector<pair<string, int> >& vFields = getIndexHeader().getFields();
    int numFields = (int)vFields.size();
    writeLine("#%%num-fields=" + Convert::toString(numFields));
    if (numFields > 0) {
        for (int i = 0; i < numFields; i++) {
            string strFieldName = vFields[i].first;
            writeLine("#%%field-" + Convert::toString(i + 1) + "=" + strFieldName);
        }
    }
    else {
        setErrMsg("Failed to write header section for index file, no fields detected.\tIndexFile: " + getIndexFile());
        bSuccessful = false;
    }
    writeLine("#");

    // write copy number states
    CallTable& callTable = getIndexHeader().getCallTable();
    writeLine("#%%max-alleles=" + Convert::toString(callTable.getMaxAlleles()));
    writeLine("#%%max-cn-states=" + Convert::toString(callTable.getMaxCNStates()));
    if (callTable.getCount() > 0) {
        // <call-name>:<call-code>:<copy-number-state>
        for (int i = 0; i < callTable.getCount(); i++) {
            string strCallName = callTable[i].getCallName();
            int iCallCode = callTable[i].getCallCode();
            int iCNState = callTable[i].getCNState();
            string strCallCodeStr = strCallName + ":" + Convert::toString(iCallCode) + ":" + Convert::toString(iCNState);
            writeLine("#%%call-code-" + Convert::toString(i + 1) + "=" + strCallCodeStr);
        }
    }
    else {
        setErrMsg("Failed to write header section for index file, no call codes detected.\tIndexFile: " + getIndexFile());
        bSuccessful = false;
    }
    writeLine("#");

    return bSuccessful;
}

// index header section
bool IndexWriter::writeHeader_v1()
{
    bool bSuccessful = true;

    vector<string>& vSampleNames = getIndexHeader().getSampleNames();
    int numSamples = (int)vSampleNames.size();

    writeLine("#%%num-probesets=" + Convert::toString(getIndexHeader().getNumProbeSets()));
    writeLine("#");
    writeLine("#%%num-samples=" + Convert::toString(numSamples));
    writeLine("#");
    writeLine("#%%format-version=" + getIndexHeader().getFormatVersionStr());
    writeLine("#%%brlmmp-transform=" + getIndexHeader().getTransformMethod());

    // write sample fields
    vector<pair<string, int> >& vFields = getIndexHeader().getFields();
    int numFields = (int)vFields.size();
    if (numFields > 0) {
        for (int i = 0; i < numFields; i++) {
            string strFieldName = vFields[i].first;
            writeLine("#%%field-" + Convert::toString(i + 1) + "=" + strFieldName);
        }
    }
    else {
        setErrMsg("Failed to write header section for index file, no fields detected.\tIndexFile: " + getIndexFile());
        bSuccessful = false;
    }

    // write sample names
    if (numSamples > 0) {
        int count = 1;
        for (int i = 0; i < numSamples; i++) {
            string strSampleName = vSampleNames[i];
            writeLine("#%%cel-name-" + Convert::toString(count) + "=" + strSampleName);
            count++;
        }
    }
    else {
        setErrMsg("Failed to write header section for index file, no samples detected.\tIndexFile: " + getIndexFile());
        bSuccessful = false;
    }
    writeLine("#");

    return bSuccessful;
}

// index body section
bool IndexWriter::writeBody()
{
    bool bSuccessful = true;

    IndexTable& indexTable = getIndexTable();
    if (indexTable.getCount() > 0) {
        indexTable.sortByRowIndex();
        writeLine("ProbeSetId\tRowIndex\tOffset" + string(getFormatVersion() >= 2.0 ? "\tNumAlleles" : ""));
        for (int i = 0; i < indexTable.getCount(); i++) {
            IndexRow& row = indexTable[i];
            if (row.getRowIndex() < 0) {
                setErrMsg("Failed to write body section for index file, invalid row index detected. Should be 0 or greater.\tRowIndex: " + Convert::toString(row.getRowIndex()) + ", IndexFile: " + getIndexFile());
                bSuccessful = false;
                break;
            }
            if (i == 0) {
                if (row.getRowIndex() != 0) {
                    setErrMsg("Failed to write body section for index file, invalid start row index. Should be 0.\tStartRowIndex: " + Convert::toString(row.getRowIndex()) + ", IndexFile: " + getIndexFile());
                    bSuccessful = false;
                    break;
                }
                if (row.getFileOffset() != 0) {
                    setErrMsg("Failed to write body section for index file, invalid start file offset. Should be 0.\tStartFileOffset: " + Convert::toString((size_t)row.getFileOffset()) + ", IndexFile: " + getIndexFile());
                    bSuccessful = false;
                    break;
                }
            }
            writeLine(row.getProbeSetId() + "\t" + Convert::toString(row.getRowIndex()) + "\t" + Convert::toString((size_t)row.getFileOffset()) + (getFormatVersion() >= 2.0 ? "\t" + Convert::toString(row.getNumAlleles()) : ""));
        }
    }
    else {
        setErrMsg("Failed to write body section for index file, no probesets detected.\tIndexFile: " + getIndexFile());
        bSuccessful = false;
    }

    return bSuccessful;
}

bool IndexWriter::initFields(int maxAlleles)
{
    try {
        clearFields();

        addField("Call:byte", sizeof(uint8_t));
        addField("Confidence:float32", sizeof(float));
        addField("LogRatio:float32", sizeof(float));
        addField("Strength:float32", sizeof(float));
        for (int i = 0; i < maxAlleles; i++) {
            string strField = "";
            strField += ('A' + i);
            strField += "AlleleSignal:float32";
            addField(strField, sizeof(float));
        }

        getIndexHeader().setMaxAlleles(maxAlleles);
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::addField(std::string field, int size)
{ 
    try {
        vector<pair<string, int> >& vFields = getIndexHeader().getFields();
        vFields.push_back(std::make_pair(field, size));
        setNumFields((int)vFields.size());
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::initSampleNames(int numSamples)
{
    try {
        clearSampleNames();
        int numFakeSamples = numSamples * getSampleMultiplier();
        getIndexHeader().getSampleNames().reserve(numFakeSamples);
        setNumSamples(0);
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::addSampleName(std::string strSampleName)
{ 
    try {
        for (int i = 0; i < getSampleMultiplier(); i++) {
            string strTempSampleName = Fs::basename(strSampleName);
            strTempSampleName += (i == 0 ? "" : "_" + Convert::toString(i));
            getIndexHeader().getSampleNames().push_back(strTempSampleName);
            incrementNumSamples();
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::initCallCodes(int numCallCodes)
{
    try {
        clearCallCodes();
        getCallTable().reserve(numCallCodes);
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::addCallCode(string strCallName, int iCallCode, int iCNState)
{
    try {
        CallRow row;
        row.setCallName(strCallName);
        row.setCallCode(iCallCode);
        row.setCNState(iCNState);
        if (!getCallTable().add(strCallName, iCallCode, iCNState)) {
            return false;
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::generateCallCodes()
{
    if (getFormatVersion() == 1.0) {
        return generateCallCodes_v1();
    }

    try {
        int maxAlleles = getMaxAlleles();
        clearCallCodes();
        addCallCode("NoCall", 0, 0);
        addCallCode("NoCall", 1, 1);
        addCallCode("NoCall", 2, 2);
        addCallCode("OTV", 3, 0);
        addCallCode("ZeroCN", 4, 0);
        int iCallCode = 5;
        for (int i = 0; i < getMaxAlleles(); i++) {
            string strCallName = "";
            strCallName += ('A' + i);
            addCallCode(strCallName, iCallCode, 1); // copy number state 1
            iCallCode++;
        }
        for (int i = 0; i < getMaxAlleles(); i++) {
            for (int j = i; j < getMaxAlleles(); j++) {
                string strCallName = "";
                strCallName += ('A' + i);
                strCallName += ('A' + j);
                addCallCode(strCallName, iCallCode, 2); // copy number state 2
                iCallCode++;
            }
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::generateCallCodes_v1()
{
    try {
        clearCallCodes();
        addCallCode("AA", 6, 0);
        addCallCode("BB", 7, 0);
        addCallCode("AB", 8, 0);
        addCallCode("NoCall", 11, 0);
        addCallCode("OTV", 12, 0);
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::initProbeSets(int numProbeSets)
{
    try {
        clearProbeSets();
        getIndexTable().reserve(numProbeSets);
        setNumProbeSets(0);
    }
    catch (...) {
        return false;
    }
    return true;
}

bool IndexWriter::addProbeSet(string strProbeSetId, int numAlleles, __int64 fileOffset, int rowIndex)
{
    try {
        IndexRow row;
        row.setProbeSetId(strProbeSetId);
        row.setFileOffset(fileOffset);
        row.setNumAlleles(numAlleles);
        row.setRowIndex(rowIndex);
        getIndexTable().push_back(row);
        incrementNumProbeSets();
    }
    catch (...) {
        return false;
    }
    return true;
}
