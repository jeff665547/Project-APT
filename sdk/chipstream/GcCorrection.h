////////////////////////////////////////////////////////////////
//
// Copyright (C) 2014 Affymetrix, Inc.
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
 * @file   GcCorrection.h
 * @author Vincent Bressler
 * @date   Feb 27, 2014
 * 
 * @brief Class for doing an signal intensity correction on all probes based
 * on the median intensity of probes with similar GC content.
 */
#ifndef _GCCR_H_
#define _GCCR_H_

#include "chipstream/PsBoard.h"
//
#include "chipstream/ChipStream.h"
#include "chipstream/ProbeSet.h"
#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Verbose.h"
//
#include <map>
#include <vector>
//

/// String describing gc adjust algorithm/module
#define GCCR "gc-correction"

/**
 *  Class for doing an adjustment based on the median intensity of probes
 * with similar GC content.
 */
class GcCorrection : public ChipStream {

public:

  /** Constructor. */
  GcCorrection(int algo_rev, bool cel_out, int maxGc=26) {
	m_algo_rev = algo_rev;
	m_cel_out = cel_out;
    m_MaxGc = maxGc;
    m_Type = GCCR;
    m_ChipCount = 0;
    setupSelfDoc(*this);
	
	m_probe_count_per_gcc = 0;
	m_probe_intensities_per_gcc = 0;
	m_ChipLayout_obj = 0;
  }

  ~GcCorrection ()
  {
	  clear_gcc_arrays ();
  }

  void compute_gcc_arrays (std::vector<float>& data);

  void clear_gcc_arrays  (void)
  {
	  if (m_probe_intensities_per_gcc)
	  {
		  for (int i=0; i<=m_MaxGc; i++)
		  {
			  if (m_probe_intensities_per_gcc[i])
			  {
				  delete [] m_probe_intensities_per_gcc[i];
				  m_probe_intensities_per_gcc[i] = 0;
			  }
		  }

		  m_probe_intensities_per_gcc = 0;
	  }

	  delete [] m_probe_count_per_gcc;
	  m_probe_count_per_gcc = 0;
  }

  /** 
   * Caluculates medians for different levels (bins) of GC count
   * between 0 and binCount in size and appends to the chipBins. The
   * probes vector contains a list of probes that are thought to be
   * representative of background at different GC counts. The
   * probeGcVec contains the GC count for all probes on the array.
   * 
   * @param data - One microarray's worth of data.
   * @param binCount - Number of bins we are calculating GC content for.
   * @param chipBins - Matrix of chip by gc count containing the
   * median intensity for a collection of gc background
   * probes. Ordering is chipBins[chip][gcCount] and parameters from
   * latest chip in data vector will be appended as last row.
   * @param probes - Probes to use for estimating background for particular GC bings.
   * @param probeGcVec - GC count of every probe on microarray indexed by id.
   * @param warnings - Print warnings if a particular bin is empty?
   * For example: print a warning if there are no probes in the GC
   * count zero bin.
   */
  static void learnParameters(const std::vector<float> &data, 
                              int binCount, 
                              std::vector<vector<float> > &chipBins,
                              /*std::vector<int> &probes, */
                              std::vector<char> &probeGcVec, 
                              bool warnings);

  /** 
   * Give the background adjusted intensity for a particular probe on a particular chip.
   * @param probeIx - Probe of interest.
   * @param chipIx - Set of Chips of interest.
   * @param intensity - Current set of intensities.
   * @param return - background adjusted set of intensities.
   */
  float transform(int probeIx, int chipIx, float intensity, PsBoard &board) {
    return transform(probeIx, chipIx, intensity);
  }

  float transform(int probeIx, int chipIx, float intensity);

  float transform(int probeIx, int chipIx, float intensity, int probeGc);

  void transform(int chipIx, std::vector<float>& intensity);

  bool is_filtered_probe (int iprobe);

  void newChip(std::vector<float> &data);

  /** 
   * @brief Method for being passed a new cel file worth of data.
   * @param data - vector of vectors of cel file data from same sample.
   */
  void newDataSet(IntensityMart* iMart);

  /** 
   * Fill in the information for Self documentation.
   * @param doc - Self documenter to be filled in.
   */
  static void setupSelfDoc(SelfDoc &doc) {
    doc.setDocName(GCCR);
    doc.setDocDescription("Correct feature intensity for variations in gc_count.");

    doc.setDocOptions(getDefaultDocOptions());
  }

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() { 
    std::vector<SelfDoc::Opt> opts;

	SelfDoc::Opt algo_rev = {"algo_rev", SelfDoc::Opt::Integer, "4", "4", "1", "4", "1 - match median gc_count=N intensities to set to median gc_count=12 intensities, 4 - match gc_count=N intensity profile to gc_count=12 intensity profile."};
	  opts.push_back(algo_rev);

	SelfDoc::Opt cel_out = {"cel_out", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", "true to create a cel file containing GcCorrected intensities."};
	opts.push_back(cel_out);

    return opts;
  }

  /** 
   * @brief Supply a little how/what/why about the algorithms this
   * class performs and what parameters it takes.
   * @return SelfDoc
   */
  static SelfDoc explainSelf() { 
    SelfDoc doc;
    setupSelfDoc(doc);
    return doc;
  }

  /** 
   * @brief This static function should be overridden by child classes
   * to return an object of the correct type initialized correctly
   * with the parameters in the string, string map. All objects
   * created this way should be deleted when finished using.
   * 
   * @param param - Map of key/value pairs to initialize the object.
   * 
   * @return Pointer toCreate object, this should be sub casted as necessary.
   */
  static SelfCreate *newObject(std::map<std::string,std::string> &param) 
  {
	SelfDoc doc = explainSelf();

	int algo_rev = 4;
	bool cel_out = false;
  
	fillInValue(algo_rev, "algo_rev", param, doc);
	fillInValue(cel_out, "cel_out", param, doc);

	GcCorrection *gccorrection = new GcCorrection(algo_rev,cel_out);
	return gccorrection;
  }

  /** 
   * Setting the GC counts for all the probes on the array.
   * @param vec - GC count for every probe on array indexed by id.
   */
  void setProbeGcVec(const std::vector<char> &vec, ChipLayout *ChipLayout_obj=0);

  /** 
   * Setting the GC counts for all the probes on the array.
   * @param vec - GC count for every probe on array indexed by id.
   */
  void setProbeGcVec(const std::vector<unsigned char> &vec);

  /** 
   * Set the probes to use for estimating background.
   * @param vec - Vector of probes that are thought to be
   * representative of background binding.
   */
  void setControlProbes(std::vector<Probe *> &vec);

  /** 
   * Set the indexes of the probes to use for estimating background.
   * @param vec - Vector of probe indexes that are thought to be
   * representative of background binding.
   */
  void setControlProbes(std::vector<int> &vec);

  /**
   * Setup the background probes and background gc indexes by reading
   * them from the blackboard
   */
  void setParameters(PsBoard &board);

private:

  /// Parameters learned (median of GC background probes).
  std::vector<vector<float> > m_Bins;
  /// Probe to be used for estimating background median
  std::vector<int> m_BgProbes;
  /// GC count of each probe matching that above
  std::vector<char> m_BgProbeGc;
  /// mdsum of probe ids
  std::string m_ProbeMd5Sum;
  /// How many chips have been seen.
  int m_ChipCount;
  /// Algorithm revision - valid values are 1 and 4
  /// algo_rev = 1 - set the median intensity for each gc count equal to the median intensity of gc count = 12
  /// algo_rev = 4 - set the intensity profile for each gc count equal to the intensity profile at gc count = 12 
  int m_algo_rev;
  /// cel_out - if true, create a cel file containing GcCorrected intensities.
  bool m_cel_out;
  /// Maximum GC count allowed.
  int m_MaxGc;

  int *m_probe_count_per_gcc;
  float **m_probe_intensities_per_gcc;

  ChipLayout *m_ChipLayout_obj;
};

#endif /* _GCCR_H_ */
