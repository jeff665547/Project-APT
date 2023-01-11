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
// ~/apt2-src/util-package/IndexTable.h ---
//
// Revision History:
// Created by David Le on 11/30/2015
//

#ifndef __INDEXTABLE_H__
#define __INDEXTABLE_H__

#include "PackageConsts.h"

#include <vector>

enum SortMode {
    BY_ROWINDEX = 1,
    BY_PROBESETID = 2
};

class IndexTable;
class IndexRow;

class IndexTable : public std::vector<IndexRow>
{
public:
    IndexTable();
    ~IndexTable();

public:
    std::string getErrMsg() { return m_strErrMsg; }
    void setErrMsg(std::string str) { m_strErrMsg = str; }
    void clearErrMsg() { m_strErrMsg = ""; }

    std::string getIndexFile() { return m_strIndexFile; }
    void setIndexFile(std::string str) { m_strIndexFile = str; }

    SortMode getSortMode() { return m_eSortMode; }
    void setSortMode(SortMode e) { m_eSortMode = e; }
    bool isSortedByRowIndex() { return (m_eSortMode == BY_ROWINDEX); }
    bool isSortedByProbeSetId() { return (m_eSortMode == BY_PROBESETID); }

    int getCount() { return (int)size(); }

public:
    void reset();
    bool readFile(std::string strIndexFile);
    bool addRow(std::string strProbeSetId, int iRowIndex, __int64 iFileOffset, int iNumAlleles = 2);
    bool isRowIndexAssigned();
    void assignRowIndexes();
    void sortByRowIndex();
    void sortByProbeSetId();
    IndexTable::iterator find(std::string strProbeSetId);

private:
    std::string m_strErrMsg;
    std::string m_strIndexFile;
    SortMode m_eSortMode;
};

class IndexRow
{
public:
    IndexRow();
    ~IndexRow();

    static bool compareRowIndex(const IndexRow& elem1, const IndexRow& elem2);
    static bool compareProbeSetId(const IndexRow& elem1, const IndexRow& elem2);

public:
    int getRowIndex() { return m_iRowIndex; }
    void setRowIndex(int i) { m_iRowIndex = i; }

    std::string getProbeSetId() { return m_strProbeSetId; }
    void setProbeSetId(std::string str) { m_strProbeSetId = str; }

    __int64 getFileOffset() { return m_iFileOffset; }
    void setFileOffset(__int64 i) { m_iFileOffset = i; }

    int getNumAlleles() { return m_iNumAlleles; }
    void setNumAlleles(int i) { m_iNumAlleles = i; }

public:
    void reset();
    std::string getRowStr();

private:
    int m_iRowIndex;
    std::string m_strProbeSetId;
    __int64 m_iFileOffset;
    int m_iNumAlleles;
};

#endif
