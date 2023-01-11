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
/**
 * @file   CallsFileReader.cpp
 * @author Anthony Kosky
 * @date   Feb 2, 2016  
 * @brief  Simple class for reading tsv file of calls generated from Axiom blob - not
 *      using TranslationInputTsvTable because of varying number of columns
 */


#include "CallsFileReader.h"
#include <stdlib.h>

const std::string CallsFileReader::s_sEmpty = ""; 

CallsFileReader::CallsFileReader() {
  m_iExp = -1;
  m_iRow = -1;
  m_iNumExp = 0;
  m_iNumRows = 0;
  m_iMaxAlleles = 2;
}

CallsFileReader::~CallsFileReader () {
  if (m_fstream.is_open()) {
    m_fstream.close();
  }
}

// readFile
//   read a Calls file generated from Axiom package using apt-package-to-calls
//
bool CallsFileReader::readFile(const std::string& sFilename) {
  m_fstream.open(sFilename.c_str(), std::fstream::in | std::fstream::binary);
  if (!m_fstream.is_open() || !m_fstream.good() || m_fstream.eof()) {
    // Need some error handling and reporting here
    return false;
  }
  // get first line
  std::string line;
  getLine(line);
  // strip any comment lines
  while (!m_fstream.eof() && (line.empty() || line[0] == '#')) {
    if (line.find("#%MaxAlleles") != std::string::npos) {
      std::string::size_type epos = line.find("=");
      std::string sMaxAlleles = line.substr(epos+1);
      m_iMaxAlleles = atoi(sMaxAlleles.c_str());
    }
    getLine(line);
  }
  if (line.empty() || line[0] == '#') {
    // Report error
    return false;
  }
  // parse header line
  std::vector<std::string> vFields;
  splitLine(line, vFields);
  int iNumFields = vFields.size();
  if (iNumFields < 1 || vFields[0] != "Probeset") {
    // Report error
    return false;
  }
  // set experiment names
  for (int i=1; i<vFields.size(); ++i) {
    m_vExperiments.push_back(vFields[i]);
  }
  m_iNumExp = m_vExperiments.size();
  m_iNumRows = 0;
  //
  // read rows
  //
  getLine(line);
  while (!(m_fstream.eof() || m_fstream.bad() || m_fstream.fail())) {
    splitLine(line, vFields);
    if (vFields.size() != iNumFields) {
      // report error
      return false;
    }
    m_vProbesetIds.push_back(vFields[0]);
    std::vector<std::string> vCalls (m_iNumExp);
    for (int i=1; i<vFields.size(); ++i) {
      vCalls.at(i-1) = vFields.at(i);
    }
    m_vvCalls.push_back(vCalls);
    m_iNumRows++;
    getLine(line);
  }
  return true;
};

bool CallsFileReader::nextRow() {
  return (++m_iRow < m_iNumRows);
}
  
bool CallsFileReader::nextExperiment() {
  m_iRow = -1; // need to call nextRow() to get first row
  return (++m_iExp < m_iNumExp);
}


bool CallsFileReader::getLine(std::string& line) {
 line.clear();
 if (!m_fstream.is_open() || !m_fstream.good() || m_fstream.eof()) {
   return false;
 }
 while (true) {
    int c = m_fstream.get();
    if (c == -1) { // eof
      break;
    }
    // CR,  LF, CR&LF, LF&CR,
    if ((c==0x0a)||(c==0x0d)) {
      // check to see if we gobble the next char too.
      int c_next=m_fstream.get();
      if (! ( ((c==0x0a)&&(c_next==0x0d)) ||
      ((c==0x0d)&&(c_next==0x0a)) )) {
        m_fstream.unget();
      }
      break;
    }
    line.push_back(c);
  }
  return true;
}

void CallsFileReader::splitLine(const std::string& line, std::vector<std::string>&  vFields) {
  vFields.clear();
  std::string::size_type pos1=0;
  std::string::size_type pos2=0;
  pos2 = line.find("\t",pos1);
  while (pos2 != std::string::npos) {
    vFields.push_back(line.substr(pos1,pos2-pos1));
    pos1 = pos2+1;
    pos2 = line.find("\t",pos1);
  }
  if (pos1 < line.size()) {
    vFields.push_back(line.substr(pos1, line.size()-pos1));
  }
}
