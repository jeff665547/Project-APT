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
// ~/apt2-src/util-package/BlobReader.cpp ---
//
// Revision History:
// Created by David Le on 12/01/2015
//

#include "BlobReader.h"
#include "util/Convert.h"
#include "util/Fs.h"
#include <iostream>
#include <sstream>
#include <string.h>
#include <memory>

#ifdef __APPLE__
#define fseeko64 fseeko
#define ftello64 ftello
#endif

using namespace std;

BlobReader::BlobReader()
{
    reset();
}

BlobReader::~BlobReader()
{
    close();
}

//
// Note: Must deallocate szData buffer after use
//
char* BlobReader::readFile(string strBlobFile, size_t& numBytes)
{
    numBytes = -1;
    size_t fileSize = BlobReader::getFileSize(strBlobFile);
    if (fileSize < 0) {
        return NULL;
    }

    FILE* fp = fopen(strBlobFile.c_str(), "rb");
    if (fp == NULL) {
        return NULL;
    }

    char* szData = new char[fileSize + 1];
    szData[0] = '\0';
    numBytes = fread(szData, 1, fileSize, fp);
    if (numBytes < 0) {
        return NULL;
    }
    szData[numBytes] = '\0';

    fclose(fp);

    return szData;
}

size_t BlobReader::getFileSize(string strFile)
{
    size_t fileSize = 0;

    FILE* file = fopen(strFile.c_str(), "rb");
    if (!file) { return 0; }

#ifdef WIN32
    _fseeki64(file, 0, SEEK_END);
    fileSize = _ftelli64(file);
#else
    fseeko64(file, 0, SEEK_END);
    fileSize = ftello64(file);
#endif

    fclose(file);

    return fileSize;
}

bool BlobReader::is64BitFile(string strFile)
{
    size_t fileSize = BlobReader::getFileSize(strFile);
    return (fileSize >= 0xffffffff);
}

int BlobReader::bytesForField(string strFieldType)
{
    if (strFieldType == "byte" || strFieldType == "char" || strFieldType == "int8") {
        return 1;
    }
    else if (strFieldType == "int16") {
        return 2;
    }
    else if (strFieldType == "int32" || strFieldType == "float32") {
        return 4;
    }
    else if (strFieldType == "int64" || strFieldType == "float64") {
        return 8;
    }
    else {
        return 0;
    }
}

void BlobReader::reset()
{
    m_strErrMsg = "";
    m_strBlobFile = "";
    m_strIndexFile = "";

    m_fileHandle = NULL;
    m_iFileSize = 0;
    m_iFilePosition = 0;
    m_bOpen = false;

    resetBuffer();
    m_indexReader.reset();
}

void BlobReader::resetBuffer()
{
    m_szBuffer[0] = '\0';
    m_iBufferLength = 0;
    m_iBufferPosition = 0;
    m_bEndOfFile = false;
}

bool BlobReader::initFolder(string strPackageFolder)
{
    if (!Fs::dirExists(strPackageFolder)) {
        setErrMsg("Failed to initialize blob reader, package folder does not exist.\tFolder: " + strPackageFolder);
        return false;
    }

    setBlobFile(Fs::join(strPackageFolder, BAT_AXIOMDATA_DIR, sBAT_BLOB_FILE));
    setIndexFile(Fs::join(strPackageFolder, BAT_AXIOMDATA_DIR, sBAT_BLOB_INDEX_FILE));

    return true;
}

bool BlobReader::init(string strBlobFile, string strIndexFile)
{
    if (strIndexFile == "") {
        strIndexFile = Fs::join(Fs::dirname(strBlobFile), sBAT_BLOB_INDEX_FILE);
    }

    if (!Fs::fileExists(strBlobFile)) {
        setErrMsg("Failed to initialize blob reader, blob file does not exist.\tBlobFile: " + strBlobFile);
        return false;
    }

    if (!Fs::fileExists(strIndexFile)) {
        setErrMsg("Failed to initialize blob reader, index file does not exist.\tIndexFile: " + strIndexFile);
        return false;
    }

    setBlobFile(strBlobFile);
    setIndexFile(strIndexFile);

    return true;
}

bool BlobReader::open()
{
    return open(getBlobFile(), getIndexFile());
}

bool BlobReader::open(string strBlobFile, string strIndexFile)
{
    clearErrMsg();

    try {
        setOpen(false);

        // load blob index file if available
        if (!readIndexFile(strIndexFile)) {
            return false;
        }

        // open blob file seeking
        setFileSize(BlobReader::getFileSize(strBlobFile));
        setFilePosition(0);
        m_fileHandle = fopen(strBlobFile.c_str(), "rb");
        if (m_fileHandle == NULL) {
            setErrMsg("Failed to open blob file for reading, file may be opened.\tBlobFile: " + strBlobFile);
            return false;
        }
    }
    catch (...) {
        setErrMsg("Exception occurred while opening blob file.\tBlobFile: " + strBlobFile);
        return false;
    }

    // blob reader is opened
    setBlobFile(strBlobFile);
    setIndexFile(strIndexFile);
    setOpen(true);

    return true;
}

bool BlobReader::close()
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

string BlobReader::readLine()
{
    setErrMsg("");

    if (!isOpen()) {
        setErrMsg("Failed to read line from stream, blob file not opened.\tBlobFile: " + getBlobFile());
        return "";
    }

    if (isEndOfFile()) {
        return "";
    }

    return nextBufferLine();
}

string BlobReader::nextBufferLine()
{
    if (m_iBufferLength == 0) {
        loadBufferData();
    }

    char szTemp[MAX_BUFFER + 1];
    string line = "";
    while (!isEndOfFile()) {
        int startPosition = getBufferPosition();
        bool endOfLine = false;
        for (; getBufferPosition() < getBufferLength(); incrementBufferPosition()) {
            if (m_szBuffer[getBufferPosition()] == '\n') {
                endOfLine = true;
                break;
            }
        }
        if (getBufferPosition() > startPosition) {
            int numBytes = getBufferPosition() - startPosition;
            memcpy(szTemp, &m_szBuffer[startPosition], numBytes);
            szTemp[numBytes] = '\0';
            line += szTemp;
        }
        if (endOfLine || isEndOfFile()) {
            if (endOfLine && line != "" && line.at(line.length() - 1) == '\r') {
#ifdef WIN32
                line.pop_back();
#else
                line = line.substr(0, line.length() - 1);
#endif
            }
            incrementBufferPosition();
            return line;
        }
        loadBufferData();
    }

    return "";
}

bool BlobReader::loadBufferData()
{
    size_t numBytes = fread(m_szBuffer, 1, MAX_BUFFER, m_fileHandle);
    if (numBytes != MAX_BUFFER) {
        setEndOfFile(true);
    }
    if (numBytes > 0) {
        setBufferLength((int)numBytes);
        setBufferPosition(0);
        return true;
    }
    return false;
}

size_t BlobReader::getCurrentFilePosition()
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to get blob current position, reader not opened.\tBlobFile: " + getBlobFile());
        return -1;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to get blob current position, file handle is invalid.\tBlobFile: " + getBlobFile());
        return -1;
    }

    size_t position = -1;

#ifdef WIN32
    position = _ftelli64(m_fileHandle);
#else
    position = ftello64(m_fileHandle);
#endif

    return position;
}

bool BlobReader::seek(size_t fileOffset)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to seek to blob offset, reader not opened.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to seek to blob offset, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    int result = 0;

#ifdef WIN32
    result = _fseeki64(m_fileHandle, fileOffset, SEEK_SET);
#else
    result = fseeko64(m_fileHandle, fileOffset, SEEK_SET);
#endif

    return (result == 0);
}

bool BlobReader::seekToFirstProbeSet()
{
    return seek(0);
}

bool BlobReader::seekToProbeSet(string strProbeSetId)
{
    clearErrMsg();

    IndexTable& indexTable = m_indexReader.getIndexTable();
    IndexTable::iterator iterFind = indexTable.find(strProbeSetId);
    if (iterFind == indexTable.end()) {
        setErrMsg("Failed to seek to blob probeset, probeset id not found.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile() + ", ProbeSetId: " + strProbeSetId);
        return false;
    }    

    return seek(iterFind->getFileOffset());
}

bool BlobReader::seekToSample(string strProbeSetId, string strSampleName)
{
    clearErrMsg();

    IndexHeader& indexHeader = m_indexReader.getIndexHeader();
    int sampleIndex = indexHeader.findSampleIndex(strSampleName);
    if (sampleIndex == -1) {
        setErrMsg("Failed to seek to blob sample, sample name not found.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile() + ", ProbeSetId: " + strProbeSetId + ", SampleName: " + strSampleName);
        return false;
    }

    IndexTable& indexTable = m_indexReader.getIndexTable();
    IndexTable::iterator iterFind = indexTable.find(strProbeSetId);
    if (iterFind == indexTable.end()) {
        setErrMsg("Failed to seek to blob sample, probeset id not found.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile() + ", ProbeSetId: " + strProbeSetId);
        return false;
    }

    size_t probeSetOffset = iterFind->getFileOffset();
    size_t sampleOffset = sampleIndex * indexHeader.getBytesPerSample(iterFind->getNumAlleles());

    return seek(probeSetOffset + sampleOffset);
}

bool BlobReader::skipToNextSample(int numAlleles, int numSamples)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to skip to next sample, reader not opened.\tBlobReader: " + getBlobFile() + ", IndexFile: " + getIndexFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to skip to next sample, file handle is invalid.\tBlobReader: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexReader.getIndexHeader();
    int bytesPerSample = numSamples * indexHeader.getBytesPerSample(numAlleles);
    int result = 0;

#ifdef WIN32
    result = _fseeki64(m_fileHandle, bytesPerSample, SEEK_CUR);
#else
    result = fseeko64(m_fileHandle, bytesPerSample, SEEK_CUR);
#endif

    return (result == 0);
}

bool BlobReader::skipToNextSampleBlock(int numAlleles)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to skip to next sample block, reader not opened.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to skip to next sample block, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexReader.getIndexHeader();
    size_t bytesPerSampleBlock = indexHeader.getBytesPerSampleBlock(numAlleles);
    int result = 0;

#ifdef WIN32
    result = _fseeki64(m_fileHandle, bytesPerSampleBlock, SEEK_CUR);
#else
    result = fseeko64(m_fileHandle, bytesPerSampleBlock, SEEK_CUR);
#endif

    return (result == 0);
}

bool BlobReader::readNextSample(SampleEx& sample, int numAlleles)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to read data for next sample, reader not opened.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to read data for next sample, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexReader.getIndexHeader();
    size_t bytesPerSample = indexHeader.getBytesPerSample(numAlleles);
    size_t currentFilePosition = getCurrentFilePosition();
    size_t fileSize = getFileSize();

    if (currentFilePosition + bytesPerSample > fileSize) {
        setErrMsg("Failed to read data for next sample, bytes exceeded stream size.\tBlobFile: " + getBlobFile() + ", RemainingBytes: " + Convert::toString(fileSize - currentFilePosition) + ", Expected: " + Convert::toString(bytesPerSample));
        return false;
    }

    char szData[MAX_BUFFER];
    size_t bytesRead = fread(szData, 1, bytesPerSample, m_fileHandle);
    if (bytesRead != bytesPerSample) {
        setErrMsg("Failed to read data for next sample, read i/o error.\tBlobFile: " + getBlobFile());
        return false;
    }
    
    if (!sample.parse(szData, numAlleles)) {
        setErrMsg("Failed to parse data for next sample, data format error.\tBlobFile: " + getBlobFile() + ", NumAlleles: " + Convert::toString(numAlleles));
        return false;
    }

    return true;
}

bool BlobReader::readNextSampleBlock(SampleList& vSamples, int numAlleles)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to read data for next sample block, reader not opened.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to read data for next sample block, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexReader.getIndexHeader();
    int bytesPerSampleBlock = indexHeader.getBytesPerSampleBlock(numAlleles);
    size_t currentFilePosition = getCurrentFilePosition();
    size_t fileSize = getFileSize();

    if (currentFilePosition + bytesPerSampleBlock > fileSize) {
        setErrMsg("Failed to read data for next sample block, bytes exceeded stream size.\tBlobFile: " + getBlobFile() + ", RemainingBytes: " + Convert::toString(fileSize - currentFilePosition) + ", Expected: " + Convert::toString(bytesPerSampleBlock));
        return false;
    }

    char szData[MAX_BUFFER];
    size_t bytesRead = fread(szData, 1, bytesPerSampleBlock, m_fileHandle);
    if (bytesRead != bytesPerSampleBlock) {
        setErrMsg("Failed to read data for next sample block, read i/o error.\tBlobFile: " + getBlobFile());
        return false;
    }

    int numSamples = indexHeader.getNumSamples();
    if (!vSamples.parseBlock(szData, numSamples, numAlleles)) {
        setErrMsg("Failed to parse data for next sample block, data format error.\tBlobFile: " + getBlobFile() + ", NumAlleles: " + Convert::toString(numAlleles));
        return false;
    }

    return true;
}

bool BlobReader::readNextWindowBlock(SampleList& vSamples, int numAlleles, int windowOffset, int windowSize)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to read data for next window block, reader not opened.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to read data for next window block, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexReader.getIndexHeader();
    int numSamples = indexHeader.getNumSamples();
    int bytesPerSample = indexHeader.getBytesPerSample(numAlleles);
    int bytesPerSampleBlock = indexHeader.getBytesPerSampleBlock(numAlleles);
    size_t currentFilePosition = getCurrentFilePosition();
    size_t fileSize = getFileSize();

    if (currentFilePosition + bytesPerSampleBlock > fileSize) {
        setErrMsg("Failed to read data for next window block, bytes exceeded stream size.\tBlobFile: " + getBlobFile() + ", RemainingBytes: " + Convert::toString(fileSize - currentFilePosition) + ", Expected: " + Convert::toString(bytesPerSampleBlock));
        return false;
    }

    char szData[MAX_BUFFER];
    size_t bytesRead = fread(szData, 1, bytesPerSampleBlock, m_fileHandle);
    if (bytesRead != bytesPerSampleBlock) {
        setErrMsg("Failed to read data for next window block, read i/o error.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (!vSamples.parseWindow(szData, numSamples, numAlleles, bytesPerSample, windowOffset, windowSize)) {
        setErrMsg("Failed to parse data for next window block, data format error.\tBlobFile: " + getBlobFile() + ", NumAlleles: " + Convert::toString(numAlleles));
        return false;
    }

    return true;
}

bool BlobReader::readIndexFile(std::string strIndexFile)
{
    clearErrMsg();

    // load blob index file if available
    if (strIndexFile != "") {
        if (!m_indexReader.readFile(strIndexFile)) {
            setErrMsg(m_indexReader.getErrMsg());
            return false;
        }
    }

    return true;
}
