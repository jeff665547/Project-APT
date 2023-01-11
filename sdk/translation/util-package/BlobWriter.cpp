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
// ~/apt2-src/util-package/BlobWriter.cpp ---
//
// Revision History:
// Created by David Le on 12/02/2015
//

#include "BlobWriter.h"
#include "util/Convert.h"
#include "util/Fs.h"

#include <iostream>
#include <string>

using namespace std;

#ifdef __APPLE__
#define fseeko64 fseeko
#define ftello64 ftello
#endif

BlobWriter::BlobWriter()
{
    reset();
    initFields(2);
}

BlobWriter::BlobWriter(int maxAlleles)
{
    reset();
    initFields(maxAlleles);
}

BlobWriter::~BlobWriter()
{
    close();
}

void BlobWriter::reset()
{
    m_strErrMsg = "";
    m_strBlobFile = "";
    m_strIndexFile = "";
    m_strBatchDataDir = BAT_AXIOMDATA_DIR;
    
    m_fileHandle = NULL;
    m_iFileSize = 0;
    m_iFilePosition = 0;
    m_eOpenMode = OPENMODE_NONE;
    m_bOpen = false;
    m_iWrittenSampleBlocks = 0;
    m_iWrittenNumAlleles = 0;

    resetBuffer();
    m_blobReader.reset();
    m_indexWriter.reset();
}

void BlobWriter::resetBuffer()
{
    m_szBuffer[0] = '\0';
    m_iBufferLength = 0;
    m_iBufferPosition = 0;
    m_bEndOfFile = false;
}

bool BlobWriter::initFolder(string strPackageFolder)
{
    clearErrMsg();

    try {
        //if (!Fs::dirExists(strPackageFolder)) {
        //    setErrMsg("Failed to initialize blob writer, package folder does not exist.\tFolder: " + strPackageFolder);
        //    return false;
        //}

        Fs::ensureWriteableDirPath(Fs::join(strPackageFolder, getBatchDataDir()));
        setBlobFile(Fs::join(strPackageFolder, getBatchDataDir(), sBAT_BLOB_FILE));
        setIndexFile(Fs::join(strPackageFolder, getBatchDataDir(), sBAT_BLOB_INDEX_FILE));
    }
    catch (...) {
        setErrMsg("Unknown exception occurred while initialize folder for blob readerFailed to initialize blob writer, package folder does not exist.\tFolder: " + strPackageFolder);
        return false;
    }

    return true;
}

bool BlobWriter::init(string strBlobFile, string strIndexFile)
{
    if (strIndexFile == "") {
        strIndexFile = Fs::join(Fs::dirname(strBlobFile), sBAT_BLOB_INDEX_FILE);
    }

    setBlobFile(strBlobFile);
    setIndexFile(strIndexFile);

    return true;
}

bool BlobWriter::open()
{
    return open(getBlobFile(), getIndexFile());
}

bool BlobWriter::open(string strBlobFile, string strIndexFile)
{
    clearErrMsg();

    try {
        setOpen(false);

        if (Fs::fileExists(strBlobFile)) {
            setOpenMode(OPENMODE_UPDATE);
            setFilePosition(0);
            m_fileHandle = fopen(strBlobFile.c_str(), "r+b"); // update blob file
            if (!m_fileHandle) {
                setErrMsg("Failed to open blob file for updating, file i/o error.\tBlobFile: " + strBlobFile);
                return false;
            }
            if (!Fs::fileExists(strIndexFile)) {
                setErrMsg("Failed to open blob writer, index file does not exist.\tIndexFile: " + strIndexFile);
                return false;
            }

            // parse index file            
            if (!m_indexWriter.open(strIndexFile)) {
                setErrMsg(m_indexWriter.getErrMsg());
                return false;
            }
        }
        else {
            setOpenMode(OPENMODE_CREATE);
            setFilePosition(0);
            m_fileHandle = fopen(strBlobFile.c_str(), "w+b"); // create blob file
            if (!m_fileHandle) {
                setErrMsg("Failed to open blob file for writing, file i/o error.\tBlobFile: " + strBlobFile);
                return false;
            }
        }
    }
    catch (...) {
        setErrMsg("Exception occurred while opening blob writer for processing.\tBlobFile: " + strBlobFile + ", IndexFile: " + strIndexFile);
        return false;
    }

    // blob file is opened
    setOpen(true);

    return true;
}

bool BlobWriter::close()
{
    bool bSuccessful = true;
    __int64 blobFileSize = 0;

    clearErrMsg();

    try {
        if (isOpen()) {
            // close blob file handle
            if (m_fileHandle != NULL) {
                fclose(m_fileHandle);
                m_fileHandle = NULL;
            }

            blobFileSize = Fs::fileSize(getBlobFile());
            if (blobFileSize != getExpectedFileSize()) {
                setErrMsg("Blob file size is not consistent with expected size, blob may be corrupted.\tFileSize: " + Convert::toString((size_t)blobFileSize) + " (expected: " + Convert::toString((size_t)getExpectedFileSize()) + ")" + ", ExpectedSampleBlocks: " + Convert::toString((size_t)getExpectedSampleBlocks()) + ", WrittenNumAlleles: " + Convert::toString((size_t)getWrittenNumAlleles()));
                bSuccessful = false;
            }
            if (getFilePosition() != getExpectedFileSize()) {
                setErrMsg("Blob file position is not consistent with expected size, blob may be corrupted.\tFilePosition: " + Convert::toString((size_t)getFilePosition()) + " (expected: " + Convert::toString((size_t)getExpectedFileSize()) + ")" + ", ExpectedSampleBlocks: " + Convert::toString((size_t)getExpectedSampleBlocks()) + ", WrittenNumAlleles: " + Convert::toString((size_t)getWrittenNumAlleles()));
                bSuccessful = false;
            }
            if (getWrittenSampleBlocks() != getExpectedSampleBlocks()) {
                setErrMsg("Incorrect # of sample blocks written to blob file, blob may be corrupted.\tWritten: " + Convert::toString((size_t)getWrittenSampleBlocks()) + " (expected: " + Convert::toString((size_t)getExpectedSampleBlocks()) + ")" + ", NumProbeSets: " + Convert::toString(getNumProbeSets()) + ", NumSamples: " + Convert::toString(getNumSamples()));
                bSuccessful = false;
            }

            // create index file
            if (!m_indexWriter.createFile(getIndexFile())) {
                setErrMsg(m_indexWriter.getErrMsg());
                bSuccessful = false;
            }

            setOpen(false);
        }
        else {
            blobFileSize = Fs::fileSize(getBlobFile());
        }

        if (bSuccessful) {
            reset();
            setFileSize(blobFileSize);
        }
    }
    catch (...) {
        setErrMsg("Exception occurred while closing blob writer.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile());
        return false;
    }

    return bSuccessful;
}

bool BlobWriter::write(const char* szData, size_t numBytes)
{
    clearErrMsg();

    try {
        if (!isOpen()) {
            setErrMsg("Failed to write data to stream, writer not opened.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile());
            return false;
        }

        if (m_fileHandle == NULL) {
            setErrMsg("Failed to write data to stream, file handle is invalid.\tBlobFile: " + getBlobFile() + ", NumBytes: " + Convert::toString(numBytes));
            return false;
        }

        size_t bytesWritten = fwrite(szData, 1, numBytes, m_fileHandle);
        if (bytesWritten != numBytes) {
            setErrMsg("Failed to write data to blob file, file i/o error.\tBlobFile: " + getBlobFile() + ", NumBytes: " + Convert::toString(numBytes));
            return false;
        }

        incrementFilePosition(numBytes);
    }
    catch (...) {
        setErrMsg("Exception occurred while writing data to stream.\tBlobFile: " + getBlobFile() + ", NumBytes: " + Convert::toString(numBytes));
        return false;
    }

    return true;
}

bool BlobWriter::writeStr(string str)
{
    return write(str.c_str(), str.size());
}

bool BlobWriter::writeLine(string line)
{
    line += "\n";
    return write(line.c_str(), line.size());
}

bool BlobWriter::writeBuffer(const char* szData, int numBytes)
{
    clearErrMsg();

    try {
        if (!isOpen()) {
            setErrMsg("Failed to write data to buffer, writer not opened.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile());
            return false;
        }

        int finalPosition = m_iBufferPosition + numBytes;
        if (finalPosition < MAX_BUFFER) {
            // transfer data to buffer
            memcpy(&m_szBuffer[m_iBufferPosition], szData, numBytes);
            incrementBufferPosition(numBytes);
            incrementFilePosition(numBytes);
        }
        else {
            // flush buffer data to blob file
            int bytesOverflow = finalPosition - MAX_BUFFER;
            int bytesToCopy = numBytes - bytesOverflow;
            memcpy(&m_szBuffer[m_iBufferPosition], szData, bytesToCopy);
            size_t bytesWritten = fwrite(m_szBuffer, 1, MAX_BUFFER, m_fileHandle);
            if (bytesWritten != MAX_BUFFER) {
                setErrMsg("Failed to flush buffer to blob file, file i/o error.\tBlobFile: " + getBlobFile());
                return false;
            }

            // reset buffer and copy overflow bytes over
            resetBuffer();
            memcpy(&m_szBuffer[m_iBufferPosition], &szData[bytesToCopy], bytesOverflow);
            incrementBufferPosition(bytesOverflow);
            incrementFilePosition(bytesOverflow);
        }
    }
    catch (...) {
        setErrMsg("Exception occurred while writing data to buffer.\tBlobFile: " + getBlobFile());
        return false;
    }

    return true;
}

bool BlobWriter::writeBufferStr(string str)
{
    return writeBuffer(str.c_str(), (int)str.size());
}

bool BlobWriter::writeBufferLine(string line)
{
    line += "\n";
    return writeBuffer(line.c_str(), (int)line.size());
}

bool BlobWriter::flushBuffer()
{
    clearErrMsg();

    try {
        if (!isOpen()) {
            setErrMsg("Failed to flush buffer to stream, writer not opened.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile());
            return false;
        }

        // flush buffer data to blob file
        if (m_iBufferPosition > 0) {
            size_t bytesWritten = fwrite(m_szBuffer, 1, m_iBufferPosition, m_fileHandle);
            if (bytesWritten != m_iBufferPosition) {
                setErrMsg("Failed to flush buffer to blob file, file i/o error.\tBlobFile: " + getBlobFile());
                return false;
            }

            // reset buffer
            resetBuffer();
        }
    }
    catch (...) {
        setErrMsg("Exception occurred while flushing buffer to blob file.\tBlobFile: " + getBlobFile());
        return false;
    }

    return true;
}

size_t BlobWriter::getCurrentFilePosition()
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to get blob current position, writer not opened.\tBlobFile: " + getBlobFile());
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

bool BlobWriter::seek(size_t fileOffset)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to seek to blob offset, writer not opened.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (!isOpenModeUpdate()) {
        setErrMsg("Failed to seek to blob offset, not in update mode.\tBlobFile: " + getBlobFile());
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

    setFilePosition(fileOffset);

    return (result == 0);
}

bool BlobWriter::seekToFirstProbeSet()
{
    return seek(0);
}

bool BlobWriter::seekToProbeSet(string strProbeSetId)
{
    clearErrMsg();

    IndexTable& indexTable = m_indexWriter.getIndexTable();
    IndexTable::iterator iterFind = indexTable.find(strProbeSetId);
    if (iterFind == indexTable.end()) {
        setErrMsg("Failed to seek to blob probeset, probeset id not found.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile() + ", ProbeSetId: " + strProbeSetId);
        return false;
    }

    return seek(iterFind->getFileOffset());
}

bool BlobWriter::seekToSample(string strProbeSetId, string strSampleName)
{
    clearErrMsg();

    IndexHeader& indexHeader = m_indexWriter.getIndexHeader();
    int sampleIndex = indexHeader.findSampleIndex(strSampleName);
    if (sampleIndex == -1) {
        setErrMsg("Failed to seek to blob sample, sample name not found.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile() + ", ProbeSetId: " + strProbeSetId + ", SampleName: " + strSampleName);
        return false;
    }

    IndexTable& indexTable = m_indexWriter.getIndexTable();
    IndexTable::iterator iterFind = indexTable.find(strProbeSetId);
    if (iterFind == indexTable.end()) {
        setErrMsg("Failed to seek to blob sample, probeset id not found.\tBlobFile: " + getBlobFile() + ", IndexFile: " + getIndexFile() + ", ProbeSetId: " + strProbeSetId);
        return false;
    }

    size_t probeSetOffset = iterFind->getFileOffset();
    size_t sampleOffset = sampleIndex * indexHeader.getBytesPerSample(iterFind->getNumAlleles());

    return seek(probeSetOffset + sampleOffset);
}

bool BlobWriter::skipToNextSample(int numAlleles, int numSamples)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to skip to next sample, writer not opened.\tBlobReader: " + getBlobFile() + ", IndexFile: " + getIndexFile());
        return false;
    }

    if (!isOpenModeUpdate()) {
        setErrMsg("Failed to skip to next sample, not in update mode.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to skip to next sample, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexWriter.getIndexHeader();
    int bytesPerSample = numSamples * indexHeader.getBytesPerSample(numAlleles);
    int result = 0;

#ifdef WIN32
    result = _fseeki64(m_fileHandle, bytesPerSample, SEEK_CUR);
#else
    result = fseeko64(m_fileHandle, bytesPerSample, SEEK_CUR);
#endif

    return (result == 0);
}

bool BlobWriter::skipToNextSampleBlock(int numAlleles)
{
    clearErrMsg();

    if (!isOpen()) {
        setErrMsg("Failed to skip to next sample block, writer not opened.\tBlobReader: " + getBlobFile() + ", IndexFile: " + getIndexFile());
        return false;
    }

    if (!isOpenModeUpdate()) {
        setErrMsg("Failed to skip to next sample block, not in update mode.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (m_fileHandle == NULL) {
        setErrMsg("Failed to skip to next sample block, file handle is invalid.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexWriter.getIndexHeader();
    size_t bytesPerSampleBlock = indexHeader.getBytesPerSampleBlock(numAlleles);
    int result = 0;

#ifdef WIN32
    result = _fseeki64(m_fileHandle, bytesPerSampleBlock, SEEK_CUR);
#else
    result = fseeko64(m_fileHandle, bytesPerSampleBlock, SEEK_CUR);
#endif

    return (result == 0);
}

bool BlobWriter::writeSample(SampleEx& sample)
{
    clearErrMsg();

    char szData[MAX_BUFFER];
    int sampleMultiplier = getSampleMultiplier();
    size_t numBytes = sample.format(szData, sampleMultiplier);
    if (numBytes == -1) {
        setErrMsg("Failed to format data for next sample, data format error.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (!write(szData, numBytes)) {
        setErrMsg("Failed to write data for next sample, file i/o error.\tBlobFile: " + getBlobFile() + ", NumBytes: " + Convert::toString(numBytes));
        return false;
    }

    incrementWrittenSampleBlocks(sampleMultiplier);
    incrementWrittenNumAlleles(sample.getNumAlleles() * sampleMultiplier);

    return true;
}

bool BlobWriter::writeSampleBlock(SampleList& vSamples)
{
    clearErrMsg();

    char szData[MAX_BUFFER];
    int sampleMultiplier = getSampleMultiplier();
    size_t numBytes = vSamples.formatBlock(szData, sampleMultiplier);
    if (numBytes == -1) {
        setErrMsg("Failed to format data for next sample block, data format error.\tBlobFile: " + getBlobFile());
        return false;
    }

    if (!write(szData, numBytes)) {
        setErrMsg("Failed to write data for next sample block, file i/o error.\tBlobFile: " + getBlobFile() + ", NumBytes: " + Convert::toString(numBytes));
        return false;
    }

    __int64 numAlleles = 0;
    for (int i = 0; i < (int)vSamples.size(); i++) {
        numAlleles += vSamples[i].getNumAlleles();
    }

    incrementWrittenSampleBlocks(vSamples.size() * sampleMultiplier);
    incrementWrittenNumAlleles(numAlleles * sampleMultiplier);

    return true;
}

bool BlobWriter::writeNextSample(SampleEx& sample)
{
    clearErrMsg();

    if (!isOpenModeUpdate()) {
        setErrMsg("Failed to write data for next sample, not in update mode.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexWriter.getIndexHeader();
    size_t bytesPerSample = indexHeader.getBytesPerSample(sample.getNumAlleles());
    size_t currentFilePosition = getCurrentFilePosition();
    size_t fileSize = getFileSize();

    if (currentFilePosition + bytesPerSample > fileSize) {
        setErrMsg("Failed to write data for next sample, bytes exceeded stream size.\tBlobFile: " + getBlobFile() + ", RemainingBytes: " + Convert::toString(fileSize - currentFilePosition) + ", Expected: " + Convert::toString(bytesPerSample));
        return false;
    }

    return writeSample(sample);
}

bool BlobWriter::writeNextSampleBlock(SampleList& vSamples)
{
    clearErrMsg();

    if (!isOpenModeUpdate()) {
        setErrMsg("Failed to write data for next sample block, not in update mode.\tBlobFile: " + getBlobFile());
        return false;
    }

    IndexHeader& indexHeader = m_indexWriter.getIndexHeader();
    int bytesPerSampleBlock = indexHeader.getBytesPerSampleBlock(vSamples);
    size_t currentFilePosition = getCurrentFilePosition();
    size_t fileSize = getFileSize();

    if (currentFilePosition + bytesPerSampleBlock > fileSize) {
        setErrMsg("Failed to write data for next sample block, bytes exceeded stream size.\tBlobFile: " + getBlobFile() + ", RemainingBytes: " + Convert::toString(fileSize - currentFilePosition) + ", Expected: " + Convert::toString(bytesPerSampleBlock));
        return false;
    }

    return writeSampleBlock(vSamples);
}
