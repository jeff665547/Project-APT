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
// ~/apt2-src/util-package/IndexReader.cpp ---
//
// Revision History:
// Created by David Le on 11/30/2015
//

#include "IndexReader.h"
#include "PackageConsts.h"
#include "util/Fs.h"

#include <iostream>
#include <sstream>
#include <string.h>
#include <memory>

using namespace std;

IndexReader::IndexReader()
{
    reset();
}

IndexReader::~IndexReader()
{
}

void IndexReader::reset()
{
    m_strErrMsg = "";
    m_strIndexFile = "";
    m_indexHeader.reset();
    m_indexTable.reset();
}

bool IndexReader::readFile(std::string strIndexFile)
{
    if (!m_indexHeader.readFile(strIndexFile)) {
        setErrMsg(m_indexHeader.getErrMsg());
        return false;
    }

    if (!m_indexTable.readFile(strIndexFile)) {
        setErrMsg(m_indexTable.getErrMsg());
        return false;
    }

    // save index file name
    setIndexFile(strIndexFile);

    return true;
}
