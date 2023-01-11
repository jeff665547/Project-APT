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
 * @file   CallsFileReader.h
 * @author Anthony Kosky
 * @date   Feb 2, 2016  
 * @brief  Simple class for reading tsv file of calls generated from Axiom blob - not
 *      using TranslationInputTsvTable because of varying number of columns
 */


#ifndef CALLS_FILE_READER_H
#define CALLS_FILE_READER_H

#include <string>
#include <vector>
#include <fstream>

class CallsFileReader
{
public:
  CallsFileReader ();
  ~CallsFileReader ();

private:
  std::fstream  m_fstream;
  std::vector<std::string>  m_vExperiments;  // CEL file names/column headers
  std::vector<std::string>  m_vProbesetIds;

  // m_vvCalls - vector of rows, each consisting of vector of calls for experiments
  std::vector<std::vector<std::string> >    m_vvCalls;

  int m_iExp;       // number of experiment
  int m_iRow;       // current row
  int m_iNumExp;   
  int m_iNumRows;   
  int m_iMaxAlleles;    // read from keyval pair

public:

  bool readFile(const std::string& sFilename);

  bool nextRow();   // return false on last row
  bool nextExperiment(); // return false on last experiment

  int  getNumExperiments() {return m_iNumExp;};
  int  getNumRows() {return m_iNumRows;};
  int  getMaxAlleles() {return m_iMaxAlleles;}

private:
  static const std::string s_sEmpty;

public:

  const std::string& getExperimentName() const {
    if (m_iExp < 0 || m_iExp >= m_iNumExp) {
      return s_sEmpty;
    }
    return m_vExperiments.at(m_iExp);
  }

  const std::string& getProbeset()  const{
    if (m_iRow < 0 || m_iRow >= m_iNumRows) {
      return s_sEmpty;
    }
    return m_vProbesetIds.at(m_iRow);
  }

  const std::string& getCall() const {
    if (m_iRow < 0 || m_iRow >= m_iNumRows || m_iExp < 0 || m_iExp >= m_iNumExp) {
      return s_sEmpty;
    }
    return m_vvCalls.at(m_iRow).at(m_iExp);
  }

private:
  bool getLine(std::string& line);

  void splitLine(const std::string& line, std::vector<std::string>&  vFields);
}; // CallsFileReader


#endif

