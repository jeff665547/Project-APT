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
// ~/apt2-src/util-package/IndexReader.h ---
//
// Revision History:
// Created by David Le on 11/30/2015
//

#ifndef __INDEXREADER_H__
#define __INDEXREADER_H__

#include "IndexHeader.h"
#include "IndexTable.h"
#include "PackageConsts.h"
#include "SampleTypesEx.h"
#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <set>
#include <vector>

class IndexReader
{
public:
    IndexReader();
    ~IndexReader();

public:
    std::string getErrMsg() { return m_strErrMsg; }
    void setErrMsg(std::string str) { m_strErrMsg = str; }
    void clearErrMsg() { m_strErrMsg = ""; }

    std::string getIndexFile() { return m_strIndexFile; }
    void setIndexFile(std::string str) { m_strIndexFile = str; }

    int getNumProbeSets() { return m_indexHeader.getNumProbeSets(); }
    int getNumSamples() { return m_indexHeader.getNumSamples(); }
    int getNumFieldsInSample() { return m_indexHeader.getNumFields(); }
    int getMaxAlleles() { return m_indexHeader.getMaxAlleles(); }

    IndexHeader& getIndexHeader() { return m_indexHeader; }
    IndexTable& getIndexTable() { return m_indexTable; }

public:
    void reset();
    bool readFile(std::string strIndexFile);

private:
    std::string m_strErrMsg;
    std::string m_strIndexFile;

    IndexHeader m_indexHeader;
    IndexTable m_indexTable;
};

#endif
