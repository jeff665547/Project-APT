////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////

/**
 * @file   MixedFileCheck.h
 * @brief  Class for testing whether two text files, containing
 * both text and numeric data, are equal, within tolerances
 * for numeric data.
 */
#ifndef MIXEDFILECHECK_H
#define MIXEDFILECHECK_H

//
#include "util/Convert.h"
#include "util/Fs.h"
#include "util/RegressionCheck.h"
//
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

/**
 * Class for testing whether two files, containing both text
 * and numeric data, are equal, within tolerances for numeric
 * data, ignoring line endings.
 */
class MixedFileCheck: public RegressionCheck
{
public:
  /**
   * Constructor.
   *
   * @param generatedFile File generated by application.
   * @param goldFile Comparison file, assumed to be correct.
   * @param eps Maximum accepted absolute difference in numeric values.
   *            i.e. if |generated-gold| >= eps then there is a difference.
   * @param skipDataLines Number of data lines to skip (after header lines auto-skipped).
   * @param allowedMismatch Maximum number of non-matching numeric values accepted.
   * @param frac Maximum accepted fractional difference in numeric values (not used by default).
   *            i.e. if |generated-gold| >= frac*max(|generated|,|gold|) then there is a difference.
   */
  MixedFileCheck (const std::string &generatedFile, const std::string &goldFile, const double eps,
    const unsigned int skipDataLines, const unsigned int allowedMismatch, const double frac = 0.0)
    : m_GeneratedFile (generatedFile), m_GoldFile (goldFile), m_Eps (eps), m_SkipDataLines (skipDataLines), 
      m_AllowedMismatch (allowedMismatch), m_MaxErrorsReport (-1), m_Frac (frac)
  { 
    m_Name=Fs::basename(generatedFile);
  }

  /** 
   * Utility function to set the max number of errors reported
   * @param max - maximum number of errors to report (set to -1 for no limit default behavior)
   * @return - void
   */
  void setMaxError(int max) {
    m_MaxErrorsReport = max;
  }

  /**
   * Check that the two files are the same, within tolerances.
   *
   * @param errorMsg Error message generated if the test fails.
   * @return bool Return true if files pass tests, else false.
   */
  bool check (std::string& errorMsg)
  {
    // Convert File Names
    m_GoldFile = Fs::convertToUncPath(m_GoldFile);
    m_GeneratedFile = Fs::convertToUncPath(m_GeneratedFile);

    // Open files.
    Verbose::out(1,"Reading in file: " + m_GoldFile);
    Verbose::out(1,"Reading in file: " + m_GeneratedFile);

    std::ifstream generatedStream;
    Fs::aptOpen(generatedStream, m_GeneratedFile);
    if (!generatedStream.is_open() && !generatedStream.good() )
    {
      errorMsg = "Unable to open generated file " + m_GeneratedFile;
      Verbose::out(1,errorMsg);
      return false;
    }
    std::ifstream goldStream;
    Fs::aptOpen(goldStream, m_GoldFile);
    if (!goldStream.is_open() && !goldStream.good() ) {
      errorMsg = "Unable to open gold file " + m_GoldFile;
      Verbose::out(1,errorMsg);
      return false;
    }
    unsigned int lineNumberGold = 0;  // absolute line number in gold file (used for reporting differences)
    unsigned int lineNumberGen = 0;   // absolute line number in generated file (used for reporting differences)
    unsigned int lineCount = 0;    // number of data lines loaded (excludes skipped header lines)
    unsigned int mismatchCount = 0;
    const char* lineEndings = "\r\n";
    std::string goldLine, generatedLine;
    Verbose::out(1,"Looking for differences.");
    while (! goldStream.eof() && ! goldStream.fail())
    {
      // Skip header lines -- for now
      // really need TsvFileCheck that replaces
      // MixedFileCheck and MatrixFileCheck
      do {
		  lineNumberGold++;
          getline (goldStream, goldLine);
      } while(goldLine.length() && ( goldLine[0] == '#'));
      do {
		lineNumberGen++;
        getline (generatedStream, generatedLine);
      } while(generatedLine.length() && (generatedLine[0] == '#'));
      if (generatedStream.eof() && ! goldStream.eof())
      {
        errorMsg = "The generated file, " + m_GeneratedFile
  	+ ", has fewer lines than the gold file, " + m_GoldFile;
        Verbose::out(1,errorMsg);
        return false;
      }
      // Skip header lines which need not be equal.
      if (++lineCount > m_SkipDataLines)
      {
        // Avoid line ending hassles.
        goldLine = goldLine.erase (goldLine.find_last_not_of (lineEndings) + 1);
        generatedLine = generatedLine.erase (generatedLine.find_last_not_of (lineEndings) + 1);

	// Skipping white space, convert each line to a series of strings.
    unsigned int columnNumber = 0;  // column number for values compared (used for reporting differences)
	std::istringstream goldStream (goldLine);
	std::istringstream generatedStream (generatedLine);
	std::string goldString, generatedString;
	while (goldStream >> goldString)
	{
      columnNumber++;
	  // Require the same number of whitespace delimited fields.
	  if (! (generatedStream >> generatedString))
	  {
        Verbose::out(1,"Unequal amount of whitespace delimited fields in files");
	    errorMsg = lineErrorMsg (generatedLine, goldLine);
        Verbose::out(1,errorMsg);
	    return false;
	  }
	  double goldDouble, generatedDouble;
	  bool goldSuccess, generatedSuccess;
	  goldDouble = Convert::toDoubleCheck (goldString.c_str(), &goldSuccess);
	  generatedDouble = Convert::toDoubleCheck (generatedString.c_str(), &generatedSuccess);

	  // If both fields are numeric, check for equality within the prescribed tolerance.
	  if (goldSuccess && generatedSuccess)
	  {
        // allowed absolute difference from fractional tolerance (zero by default)
        double eps2 = m_Frac*Max( fabs(goldDouble), fabs(generatedDouble) );
	    // absolute difference is acceptable if it satisfies either (least restrictive) tolerance
	    if (fabs (goldDouble - generatedDouble) > Max(m_Eps,eps2)) {
	      ++mismatchCount;
		  // report differences with ZERO-based line and column numbers
	      if( (int)mismatchCount<=m_MaxErrorsReport || m_MaxErrorsReport<0 )
            Verbose::out(1,"Numbers differ at gold line " + ToStr(lineNumberGold-1) 
			  + " generated line " + ToStr(lineNumberGen-1) 
			  + " column " + ToStr(columnNumber-1) 
			  + ": " + goldString + " and " + generatedString);
	      if( (int)mismatchCount==m_MaxErrorsReport+1 && m_MaxErrorsReport>0 )
            Verbose::out(1,"Number of differences exceeds maximum number (" + ToStr(m_MaxErrorsReport) 
			  + ") to report.");
		}
	    continue;
	  }

	  // If neither field is numeric, require them to be identical.
	  if ((! goldSuccess) && (! generatedSuccess) && (goldString == generatedString))
	    continue;
      else
        Verbose::out(1,"Strings differ: " + goldString + " and " + generatedString);

	  // Quit if there is a type mismatch or both fields are non-numeric, not identical.
	  errorMsg = lineErrorMsg (generatedLine, goldLine);
      Verbose::out(1,errorMsg);
	  return false;
	} // end while (goldStream >> goldString)

	// Require that the two lines have the same number of fields.
	if (generatedStream >> generatedString)
	{
	  errorMsg = lineErrorMsg (generatedLine, goldLine);
      Verbose::out(1,errorMsg);
	  return false;
	}
      } // end if (++lineCount > m_SkipDataLines)
    }   // end while (! goldStream.eof() && ! goldStream.fail())

    // The two files should reach eof at the same time.
    if (! generatedStream.eof())
    {
      errorMsg = "The generated file, " + m_GeneratedFile
        + ", has more lines than the gold file, " + m_GoldFile;
      Verbose::out(1,errorMsg);
      return false;
    }

    // Require that the number of numeric differences above tolerance is below
    // the defined threshold.
    if (mismatchCount > m_AllowedMismatch)
    {
      errorMsg = "There were " + ToStr (mismatchCount) + " instances where "
        + "numeric fields differed by more than the accepted tolerance: only "
        + ToStr (m_AllowedMismatch) + " are allowed";
      Verbose::out(1,errorMsg);
      return false;
    }

    Verbose::out(1,"Same.");
    return true;
  }

private:
  /**
   * Generate a generic error message for a line mismatch.
   *
   * @param generatedLine Line generated by the application.
   * @param goldLine Line considered to be correct.
   * @return Error message.
   */
  const std::string lineErrorMsg (const std::string& generatedLine, const std::string& goldLine)
  {
    const std::string msg = "Mismatch reading generated file " + m_GeneratedFile
      + ":\ngold line: '" + goldLine + "'\ngenerated line: '" + generatedLine + "'";
    return msg;
  }

  /// Name of file generated by application being tested.
  std::string m_GeneratedFile;
  /// Name of file assumed to be correct.
  std::string m_GoldFile;
  /// Maximum accepted absolute difference in numeric values.
  const double m_Eps;
  /// Number of data lines to skip (header lines are already auto-skipped).
  const unsigned int m_SkipDataLines;
  /// Maximum number of non-matching numeric values accepted.
  const unsigned int m_AllowedMismatch;
  /// The maximum number of errors to report
  int m_MaxErrorsReport;
  /// Maximum fractional difference considered equivalent.
  const double m_Frac;
};

#endif /* MIXEDFILECHECK_H */
