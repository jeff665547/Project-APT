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
// ~/apt2-src/util-package/IndexHeader.h ---
//
// Revision History:
// Created by David Le on 11/30/2015
//

#ifndef __INDEXHEADER_H__
#define __INDEXHEADER_H__

#include "CallTable.h"
#include "PackageConsts.h"
#include "SampleTypesEx.h"

#include <stdio.h>
#include <map>
#include <vector>

class IndexHeader
{
public:
    static int bytesForField(std::string strFieldType);

public:
    IndexHeader();
    ~IndexHeader();

public:
    std::string getErrMsg() { return m_strErrMsg; }
    void setErrMsg(std::string str) { m_strErrMsg = str; }
    void clearErrMsg() { m_strErrMsg = ""; }

    std::string getIndexFile() { return m_strIndexFile; }
    void setIndexFile(std::string str) { m_strIndexFile = str; }

    double getFormatVersion() { return m_dFormatVersion; }
    void setFormatVersion(double d) { m_dFormatVersion = d; }
    std::string getFormatVersionStr() { char sz[MAX_LINE]; sprintf(sz, "%.1lf", m_dFormatVersion); return sz; }

    std::string getTransformMethod() { return m_strTransformMethod; }
    void setTransformMethod(std::string str) { m_strTransformMethod = str; }

    int getNumProbeSets() { return m_iNumProbeSets; }
    void setNumProbeSets(int i) { m_iNumProbeSets = i; }
    void incrementNumProbeSets() { m_iNumProbeSets++; }

    int getNumSamples() { return m_iNumSamples; }
    void setNumSamples(int i) { m_iNumSamples = i; }
    void incrementNumSamples() { m_iNumSamples++; }

    int getNumFields() { return m_iNumFields; }
    void setNumFields(int i) { m_iNumFields = i; }

    int getFixedFields() { return m_iFixedFields; }
    void setFixedFields(int i) { m_iFixedFields = i; }

    int getSignalFields() { return m_iSignalFields; }
    void setSignalFields(int i) { m_iSignalFields = i; }

    std::vector<std::string>& getSampleNames() { return m_vSampleNames; }
    std::vector<std::pair<std::string, int> >& getFields() { return m_vFields; }
    std::map<std::string, int>& getMapSampleToIndexes() { return m_mapSampleToIndexes; }
    std::map<std::string, int>& getMapFieldToIndexes() { return m_mapFieldToIndexes; }

    int getMaxAlleles() { return m_callTable.getMaxAlleles(); }
    void setMaxAlleles(int i) { return m_callTable.setMaxAlleles(i); }

    int getMaxCNStates() { return m_callTable.getMaxCNStates(); }
    void setMaxCNStates(int i) { return m_callTable.setMaxCNStates(i); }

    int getNumCallRows() { return m_callTable.getCount(); }

    CallTable& getCallTable() { return m_callTable; }

public:
    void reset();
    bool readFile(std::string strIndexFile);
    int getBytesPerSample(int numAlleles);
    int getBytesPerSampleBlock(int numAlleles);
    int getBytesPerSampleBlock(SampleList& vSamples);
    int findSampleIndex(std::string strSampleName);

private:
    std::string m_strErrMsg;
    std::string m_strIndexFile;

    double m_dFormatVersion;
    //
    // http://www.affymetrix.com/estore/support/developer/powertools/changelog/apt-probeset-genotype.html.affx#transformations
    //   MVA = Minus vs Average (default)
    //   CES = Contrast Extremes Stretch
    //   CCS = Contrast Centers Stretch
    //   RVT = R vs Theta
    //
    std::string m_strTransformMethod;
    int m_iNumProbeSets;
    int m_iNumSamples;
    int m_iNumFields;
    int m_iFixedFields;
    int m_iSignalFields;

    std::vector<std::string> m_vSampleNames;                // <cel_name> for each sample
    std::vector<std::pair<std::string, int> > m_vFields;    // <field_name, num_bytes> -- Call:byte|Confidence:float32|LogRatio:float32|Strength:float32|AAlleleSignal:float32|BAlleleSignal:float32|CAlleleSignal:float32|... (up to 16 alleles)
    std::map<std::string, int> m_mapSampleToIndexes;        // <cel_name_lower, index>
    std::map<std::string, int> m_mapFieldToIndexes;         // <field_name_lower, index> -- call|confidence|logratio|strength

    CallTable m_callTable;
};

#endif
