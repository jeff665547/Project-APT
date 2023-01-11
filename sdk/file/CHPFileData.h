////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#ifndef _CHPFileData_HEADER_
#define _CHPFileData_HEADER_

/*! \file CHPFileData.h This file provides CHP file reading capabilities.
 */

//////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
#pragma warning(disable: 4786) // identifier was truncated in the debug information
#endif

//////////////////////////////////////////////////////////////////////

#include "file/TagValuePair.h"
//
#include "portability/affy-base-types.h"
//
#include <cstring>
#include <fstream>
#include <list>
#include <string>
#include <vector>
//

//////////////////////////////////////////////////////////////////////

namespace affxchp
{

//////////////////////////////////////////////////////////////////////

/*! This class stores a zone's background value */
typedef struct _BackgroundZoneType
{
	/*! The X coordinate of the center of the zone. */
	float centerx;

	/*! The Y coordinate of the center of the zone. */
	float centery;

	/*! The zone's background value */
	float background;

	/*! Assignment operator
	 * @param zn The zone to copy
	 * @return The new zone object
	 */
	_BackgroundZoneType operator=(_BackgroundZoneType zn)
	{
		centerx = zn.centerx; 
		centery = zn.centery; 
		background = zn.background; 
		return *this; 
	}
} BackgroundZoneType;

/*! The size of the zone as stored in a CHP file. */
#define ZONE_INFO_TYPE_SIZE (3*sizeof(float))

/*! An STL list of zones */
typedef std::list<BackgroundZoneType> BackgroundZoneTypeList;

/*! Stores a list of zones for the entire array */
typedef struct _BackgroundZoneInfo
{
	/*! The number of zones used in the array */
	int number_zones;

	/*! The smoothing factor used to calculate the zone backgrounds */
	float smooth_factor;

	/*! The list of zone background values */
	BackgroundZoneTypeList zones;
} BackgroundZoneInfo;

////////////////////////////////////////////////////////////////////

/*! This class provides storage for the CHP file header */
class CCHPFileHeader
{
public:
	/*! Constructor */
	CCHPFileHeader();

	/*! Destructor */
	~CCHPFileHeader();

public:
	/*! Defines the assay type for the array */
	typedef enum {Expression, Genotyping, Resequencing, Universal, Unknown} GeneChipAssayType;

protected:
	/*! The magic number in the file */
	int m_Magic;

	/*! The version number in the file */
	int m_Version;

	/*! The number of feature columns in the array */
	unsigned short m_Cols;

	/*! The number of feature rows in the array */
	unsigned short m_Rows;

	/*! The number of probe set results */
	int m_NumProbeSets;

	/*! The type of results stored in the CHP file */
	GeneChipAssayType m_AssayType;

	/*! The chip type or probe array type of the CHP file */
	std::string m_ChipType;

	/*! The name of the algorithm used to create the CHP file */
	std::string m_AlgorithmName;

	/*! The version number of the algorithm used to create the CHP file */
	std::string m_AlgorithmVersion;

	/*! The name of the CEL file used in the creation of the CHP file */
	std::string m_ParentCellFile;

	/*! The programmatic identifier of the algorithm used to create the CHP file */
	std::string m_ProgID;

	/*! The list of algorithm parameters */
	TagValuePairTypeList m_AlgorithmParameters;

	/*! A list of summary parameters generated by the CHP file generating algorithm */
	TagValuePairTypeList m_SummaryParameters;

	/*! The background's for each of the zones (calculated by the expression algorithm) */
	BackgroundZoneInfo m_BackgroundZoneInfo;

	/*! Parses the parameters string into a list given the delimiters.
	 * @param tagList The resulting parameter name/value list.
	 * @param strSource The parameters in a string representation
	 * @param sDelimiter1 The delimiter between each parameter
	 * @param sDelimiter2 The delimiter between the tag and value
	 */
	void ParseString(TagValuePairTypeList &tagList,
					 std::string strSource,
					 std::string sDelimiter1, 
					 std::string sDelimiter2);

	/*! Clears the class members */
	void Clear();

	/*! Friend to the parent object */
	friend class CCHPFileData;

public:
	/*! Gets the number of feature columns
	 * @return The number of feature columns
	 */
	int GetCols() const { return m_Cols; }

	/*! Gets the number of feature rows
	 * @return The number of feature rows
	 */
	int GetRows() const { return m_Rows; }

	/*! Gets the number of probe sets
	 * @return The number of probe sets
	 */
	int GetNumProbeSets() const { return m_NumProbeSets; }

	/*! Gets the assay type
	 * @return The assay type
	 */
	GeneChipAssayType GetAssayType() const { return m_AssayType; }

	/*! Gets the chip type
	 * @return The chip type
	 */
	std::string GetChipType() const { return m_ChipType; }

	/*! Gets the algorithm name
	 * @return The algorithm name
	 */
	std::string GetAlgName() const { return m_AlgorithmName; }

	/*! Gets the algorithm version
	 * @return The algorithm version
	 */
	std::string GetAlgVersion() const { return m_AlgorithmVersion; }

	/*! Gets the algorithm parameters
	 * @return The number of feature columns
	 */
	TagValuePairTypeList &AlgorithmParameters() { return m_AlgorithmParameters; }

	/*! Gets the summary parameters
	 * @return The summary parameters
	 */
	TagValuePairTypeList &SummaryParameters() { return m_SummaryParameters; }

	/*! Gets the parent CEL file
	 * @return The parent CEL file
	 */
	std::string GetParentCellFile() const { return m_ParentCellFile; }

	/*! Gets the prog ID
	 * @return The prog ID
	 */
	std::string GetProgID() const { return m_ProgID; }

	/*! Gets a specific algorithm parameter given a name/tag
	 * @return The specific algorithm parameter given a name/tag
	 */
	std::string GetAlgorithmParameter(const char *tag);

	/*! Gets a specific summary parameter given a name/tag
	 * @return The specific summary parameter given a name/tag
	 */
	std::string GetSummaryParameter(const char *tag);

	/*! Gets the background zone information
	 * @return The background zone information
	 */
	BackgroundZoneInfo &GetBackgroundZoneInfo() { return m_BackgroundZoneInfo; }

	/*! Gets the list of background zone positions and values
	 * @return The list of background zone positions and values
	 */
	BackgroundZoneTypeList &GetBackgroundZones() { return m_BackgroundZoneInfo.zones; }

	/*! Gets the background value for a given center coordinate
	 * @return The background value for a given center coordinate
	 */
	BackgroundZoneType GetBackgroundZone(int x, int y);

	/*! Gets the magic number
	 * @return The magic number
	 */
	int GetMagicNumber() const { return m_Magic; }

	/*! Gets the version number
	 * @return The version number
	 */
	int GetVersionNumber() const { return m_Version; }

	/*! Sets the number of columns
	 * @param n The number of columns.
	 */
	void SetCols(int n) { m_Cols=n; }

	/*! Sets the number of rows
	 * @param n The number of rows.
	 */
	void SetRows(int n) { m_Rows=n; }

	/*! Sets the number of probe sets
	 * @param n The number of probe sets.
	 */
	void SetNumProbeSets(int n) { m_NumProbeSets=n; }

	/*! Sets the assay type
	 * @param t The assay type.
	 */
	void SetAssayType(GeneChipAssayType t) { m_AssayType=t; }

	/*! Sets the chip type or probe array type
	 * @param s The chip type or probe array type.
	 */
	void SetChipType(const char *s) { m_ChipType=s; }

	/*! Sets the algorithm name
	 * @param s The algorithm name.
	 */
	void SetAlgName(const char *s) { m_AlgorithmName=s; }

	/*! Sets the algorithm version
	 * @param s The algorithm version.
	 */
	void SetAlgVersion(const char *s) { m_AlgorithmVersion=s; }

	/*! Sets the parent CEL file
	 * @param s The parent CEL file.
	 */
	void SetParentCellFile(const char *s) { m_ParentCellFile=s; }

	/*! Sets the prog ID
	 * @param s The programatic identifier
	 */
	void SetProgID(const char *s) { m_ProgID=s; }
};

////////////////////////////////////////////////////////////////////

/*! Provides a base class for probe set results */
class CProbeSetResults
{
public:
	/*! Constructor */
	CProbeSetResults() {};

	/*! Destructor */
	virtual ~CProbeSetResults() {};
};

////////////////////////////////////////////////////////////////////

/*! Present call for expression analysis */
#define ABS_PRESENT_CALL 0

/*! Marginal call for expression analysis */
#define ABS_MARGINAL_CALL 1

/*! Absent call for expression analysis */
#define ABS_ABSENT_CALL 2

/*! No call call for expression analysis */
#define ABS_NO_CALL 3

/*! Increase call for expression comparison analysis */
#define COMP_INCREASE_CALL 1

/*! Decrease call for expression comparison analysis */
#define COMP_DECREASE_CALL 2

/*! Moderate increase call for expression comparison analysis */
#define COMP_MOD_INCREASE_CALL 3

/*! Moderate decrease call for expression comparison analysis */
#define COMP_MOD_DECREASE_CALL 4

/*! No change call for expression comparison analysis */
#define COMP_NO_CHANGE_CALL 5

/*! No call call for expression comparison analysis */
#define COMP_NO_CALL 6

/*! Expression analysis probe set results for the MAS5 algorithm */
class CExpressionProbeSetResults : public CProbeSetResults
{
public:
	/*! The detection p-value */
	float DetectionPValue;

	/*! The signal value */
	float Signal; 

	/*! The number of probe pairs in the set */
	unsigned short NumPairs;

	/*! The number of probe pairs used to calculate the signal value */
	unsigned short NumUsedPairs;

	/*! The detection call */
	unsigned char Detection;

	/*! Flag indicating that comparison results exist */
	bool m_HasCompResults;

	/*! The change p-value */
	float ChangePValue; 

	/*! The signal log ratio */
	float SignalLogRatio; 

	/*! The signal log ratio low value */
	float SignalLogRatioLow; 

	/*! The signal log ratio high value */
	float SignalLogRatioHigh; 

	/*! The number of probe pairs in common between the experiment and baseline data */
	unsigned short NumCommonPairs;

	/*! The change call */
	unsigned char Change;

	/*! Returns a string representation of the detection call.
	 * @return The detection call
	 */
	std::string GetDetectionString();

	/*! Returns a string representation of the change call.
	 * @return The change call
	 */
	std::string GetChangeString();

	/*! Assignment operator
	 * @param src The object to copy
	 * @return The copied object
	 */
	CExpressionProbeSetResults operator=(CExpressionProbeSetResults &src);

	/*! Constructor */
	CExpressionProbeSetResults() { m_HasCompResults = false; }

	/*! Destructor */
	~CExpressionProbeSetResults() {}
};

////////////////////////////////////////////////////////////////////

/*! The AA allele call */
#define ALLELE_A_CALL 6

/*! The BB allele call */
#define ALLELE_B_CALL 7

/*! The AB allele call */
#define ALLELE_AB_CALL 8

/*! The no call allele call */
#define ALLELE_NO_CALL 11

/*! Genotyping analysis probe set results */
class CGenotypeProbeSetResults : public CProbeSetResults
{
public:
	/*! The allele call */
	unsigned char AlleleCall;

	/*! The confidence associated with the allele call */
	float Confidence;

	/*! The relative allele strength for the first block */
	float RAS1;

	/*! The relative allele strength for the second block */
	float RAS2;

	/*! The p-value associated with an AA call */
	float pvalue_AA;

	/*! The p-value associated with an AB call */
	float pvalue_AB;

	/*! The p-value associated with an BB call */
	float pvalue_BB;

	/*! The p-value associated with an no call call */
	float pvalue_NoCall;

	/*! Returns a string representation of the allele call.
	 * @return The allele call
	 */
	std::string GetAlleleCallString();

	/*! Assignment operator
	 * @param src The object to copy
	 * @return The copied object
	 */
	CGenotypeProbeSetResults operator=(CGenotypeProbeSetResults &src);

	/*! Constructor */
	CGenotypeProbeSetResults() {Confidence=RAS1=RAS2=pvalue_AA=pvalue_AB=pvalue_BB=pvalue_NoCall=0;}

	/*! Destructor */
	~CGenotypeProbeSetResults() {}
};

////////////////////////////////////////////////////////////////////

/*! Universal (tag array) analysis probe set results. */
class CUniversalProbeSetResults : public CProbeSetResults
{
protected:
	/*! The background value.*/
	float background;

public:
	/*! Gets the background value.
	 * @return The background value.
	 */
	float GetBackground() const { return background; }

	/*! Sets the background value.
	 * @param bg The background value.
	 */
	void SetBackground(float bg) { background = bg; }

	/*! Assignment operator
	 * @param src The object to copy
	 * @return The copied object
	 */
	CUniversalProbeSetResults operator=(CUniversalProbeSetResults &src);

	/*! Constructor */
	CUniversalProbeSetResults() { background=0; }

	/*! Destructor */
	~CUniversalProbeSetResults() {}
};

////////////////////////////////////////////////////////////////////

/*! A structure to hold a force call, its position and reason.
 *
 * A force call is the call the algorithm would have made if the thresholds
 * were not applied.
 */
typedef struct _ForceCallType
{
	/*! The position (index) of the call. */
	int position;

	/*! The force call. */
	char call;

	/*! The reason for the call. */
	unsigned char reason;

} ForceCallType;

/*! The force call was made due to no signal threshold. */
#define NO_SIGNAL_THR_FORCE_CALL 'N'

/*! The force call was made due to weak signal threshold. */
#define WEAK_SIGNAL_THR_FORCE_CALL 'W'

/*! The force call was made due to saturation level. */
#define SATURATION_LEVEL_FORCE_CALL 'S'

/*! The force call was made due to quality score threshold. */
#define QUALITY_SCORE_THR_FORCE_CALL 'Q'

/*! The force call was made due to failed both trace and sequence profiles. */
#define TRACE_AND_SEQUENCE_PROFILES_FORCE_CALL 'F'

/*! The force call was made due to base reliability threshold. */
#define RELIABILITY_THR_FORCE_CALL 'B'

/*! A structure to hold a base call at a given position (index). */
typedef struct _BaseCallType
{
	/*! The position (index) of the call. */
	int position;

	/*! The call. */
	char call;

} BaseCallType;

/*! Resequencing results. */
class CResequencingResults
{
protected:
	/*! The called bases. */
	std::vector<char> calledBases;

	/*! Base call scores. */
	std::vector<float> scores;

	/*! An array of force calls - base calls the algorithm would have made if the thresholds were removed. */
	std::vector<ForceCallType> forceCalls;

	/*! An array of original calls. The calledBases contained the results of the algorithm and user edits.
	 * If a user edits a base the original algorithm called base is stored in this vector.
	 */
	std::vector<BaseCallType> origCalls;

public:
	/*! Constructor */
	CResequencingResults() {}

	/*! Destructor */
	~CResequencingResults() { Clear(); }

	/*! Clears the members. */
	void Clear() { calledBases.clear(); scores.clear(); forceCalls.clear(); origCalls.clear(); }

	/*! Gets the called bases.
	 * @return The array of called bases.
	 */
	const std::vector<char> &GetCalledBases() { return calledBases; }

	/*! Gets the called base at the given index.
	 * @param index The index to the called bases array.
	 * @return The called base.
	 */
	char GetCalledBase(int index) { return calledBases[index]; }

	/*! Gets the size of the called bases array.
	 * @return The size of the called bases array.
	 */
	int GetCalledBasesSize() const { return (int) calledBases.size(); }

	/*! Resizes the called bases array.
	 * @param size The size of the array.
	 */
	void ResizeCalledBases(int size) { calledBases.resize(size); }

	/*! Sets the called base.
	 * @param index The index to the array.
	 * @param call The call.
	 */
	void SetCalledBase(int index, char call) { calledBases[index] = call; }

	/*! Gets the scores.
	 * @return The array of scores.
	 */
	const std::vector<float> &GetScores() { return scores; }

	/*! Gets the score at the given index.
	 * @param index The index to the scores array.
	 * @return The score.
	 */
	float GetScore(int index) { return scores[index]; }

	/*! Gets the size of the scores array.
	 * @return The size of the scores array.
	 */
	int GetScoresSize() const { return (int) scores.size(); }

	/*! Resizes the scores array.
	 * @param size The size of the array.
	 */
	void ResizeScores(int size) { scores.resize(size); }

	/*! Sets the score.
	 * @param index The index to the array.
	 * @param score The score.
	 */
	void SetScore(int index, float score) { scores[index] = score; }

	/*! Gets the force calls.
	 * @return The array of force calls.
	 */
	const std::vector<ForceCallType> &GetForceCalls() { return forceCalls; }

	/*! Gets the force call at the given index.
	 * @param index The index to the force calls array.
	 * @return The force call.
	 */
	ForceCallType GetForceCall(int index) { return forceCalls[index]; }

	/*! Gets the size of the force calls array.
	 * @return The size of the force calls array.
	 */
	int GetForceCallsSize() const { return (int) forceCalls.size(); }

	/*! Resizes the force calls array.
	 * @param size The size of the array.
	 */
	void ResizeForceCalls(int size) { forceCalls.resize(size); }

	/*! Sets the force call.
	 * @param index The index to the array.
	 * @param call The force call.
	 */
	void SetForceCall(int index, ForceCallType call) { forceCalls[index] = call; }

	/*! Gets the original called bases.
	 * @return The array of original calls.
	 */
	const std::vector<BaseCallType> &GetOrigCalls() { return origCalls; }

	/*! Gets the original called base at the given index.
	 * @param index The index to the original calls array.
	 * @return The original call.
	 */
	BaseCallType GetOrigCall(int index) { return origCalls[index]; }

	/*! Gets the size of the original calls array.
	 * @return The size of the original calls array.
	 */
	int GetOrigCallsSize() const { return (int) origCalls.size(); }

	/*! Resizes the original calls array.
	 * @param size The size of the array.
	 */
	void ResizeOrigCalls(int size) { origCalls.resize(size); }

	/*! Sets the original call.
	 * @param index The index to the array.
	 * @param call The original call.
	 */
	void SetOrigCall(int index, BaseCallType call) { origCalls[index] = call; }
};

////////////////////////////////////////////////////////////////////

/*! This class provides storage and reading capabilities for CHP files */
class CCHPFileData  
{
public:
	/*! Constructor */
	CCHPFileData();

	/*! Destructor */
	~CCHPFileData();

protected:
	/*! The file header object */
	CCHPFileHeader m_Header;

	/*! The full path of the CHP file */
	std::string m_FileName;

	/*! A string to hold an error message associated with a read operation */
	std::string m_strError;

	/*! The vector of probe set results */
	std::vector<CProbeSetResults *> m_ProbeSetResults;

	/*! The resequencing results. */
	CResequencingResults m_ReseqResults;

	/*! Opens the file for reading.
	 * @param bReadHeaderOnly Flag to indicate if the header is to be read only.
	 * @return True if successful.
	 */
	bool Open(bool bReadHeaderOnly=false);

public:
	/*! Accessors to header.
	 * @return The header data object
	 */
	CCHPFileHeader &GetHeader() { return m_Header; }

	/*! Returns the expression probe set result
	 * @param index The index to the result object of interest.
	 * @return The expression result.
	 */ 
	CExpressionProbeSetResults *GetExpressionResults(int index);

	/*! Returns the genotyping probe set result
	 * @param index The index to the result object of interest.
	 * @return The genotyping result.
	 */
	CGenotypeProbeSetResults *GetGenotypingResults(int index);

	/*! Returns the universal (tag array) probe set result
	 * @param index The index to the result object of interest.
	 * @return The universal result.
	 */
	CUniversalProbeSetResults *GetUniversalResults(int index);

	/*! Returns the resequencing results.
	 * @return The resequencing results.
	 */
	CResequencingResults *GetResequencingResults();

	/*! Error string when the read functions fail.
	 * @return A string message describing a read error
	 */
	std::string GetError() const { return m_strError; }

	// Functions to read file.
	bool Read();

	/*! Reads the header of the CHP file
	 * @return True if successful
	 */
	bool ReadHeader();

	/*! Determines if the file specified by the FileName property exists.
	 * @return True if the file exists.
	 */
	bool Exists();

	/*! Determines if the CHP file specified by the FileName property is an XDA format file
	 * @return True if the file is an XDA file
	 */
	bool IsXDACompatibleFile();

	bool IsMas5File();

	/*! Sets the file name.
	 * @param name The full path to the CHP file
	 */
	void SetFileName(const char *name) { m_FileName = name; }

	/*! Gets the file name.
	 * @return The full path to the CHP file.
	 */
	std::string GetFileName() const { return m_FileName; }

	/*! Deallocates any memory used by the class object */
	void Clear();
};

////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////

#endif // !defined(_CHPFileData_HEADER_)
