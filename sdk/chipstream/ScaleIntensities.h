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
 * @file   ScaleIntensities.h
 * @author Vincent Bressler
 * @date   May 20, 2014
 * 
 * @brief Scale cel file intensities
 */
#ifndef _SCIR_H_
#define _SCIR_H_

#include "chipstream/PsBoard.h"
//
#include "chipstream/ChipStream.h"
//#include "chipstream/ProbeSet.h"
//#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Verbose.h"
//
#include <map>
#include <vector>
//

/// String describing gc adjust algorithm/module
#define SCIR "scale-intensities"

/**
 *  Class for scaling cel intensities.
 */
class ScaleIntensities : public ChipStream {

public:

	ScaleIntensities(int floor, int low, int high, int ceiling, double low_pct, double high_pct, double gain_lock_pct, bool cel_out) 
	{
		m_Type = SCIR;

		m_floor = floor;
		m_low = low;
		m_high = high;
		m_ceiling = ceiling;
		m_low_pct = low_pct;
		m_high_pct = high_pct;
		m_gain_lock_pct = gain_lock_pct;

		m_cel_out = cel_out;

		setupSelfDoc(*this);

		m_ChipLayout_obj = 0;
	}

	~ScaleIntensities ()
	{
	}

	///////////////////////////////////////////////////////////////////////////////////
	// START Members that are requried by ChipStream
	static void setupSelfDoc(SelfDoc &doc) 
	{
		doc.setDocName(SCIR);
		doc.setDocDescription("Scale cel intensities.");

		doc.setDocOptions(getDefaultDocOptions());
	}

	static SelfDoc explainSelf() 
	{ 
		SelfDoc doc;
		setupSelfDoc(doc);
		return doc;
	}

	static SelfCreate *newObject(std::map<std::string,std::string> &param) 
	{
		SelfDoc doc = explainSelf();

		int floor = 1;
		int low = 20;
		int high = 50000;
		int ceiling = 1000000;

		double low_pct = 0.02;
		double high_pct = 0.98;

		double gain_lock_pct = 0;

		bool cel_out = false;

		fillInValue(floor, "floor", param, doc);
		fillInValue(low, "low", param, doc);
		fillInValue(high, "high", param, doc);
		fillInValue(ceiling, "ceiling", param, doc);

		fillInValue(low_pct, "low_pct", param, doc);
		fillInValue(high_pct, "high_pct", param, doc);

		fillInValue(gain_lock_pct,"gain_lock_pct", param, doc);

		fillInValue(cel_out, "cel_out", param, doc);

		ScaleIntensities *obj = new ScaleIntensities(floor,low,high,ceiling,low_pct,high_pct,gain_lock_pct,cel_out);
		return obj;
	}

	// iMart - contains intensity data for all the cel files in the current run
	void newDataSet(IntensityMart* iMart);

	inline void setChipLayout (ChipLayout *obj)
	{
		m_ChipLayout_obj = obj;
	}

	// END Members that are requried by ChipStream
	///////////////////////////////////////////////////////////////////////////////////

	
	void transform(int chipIx, std::vector<float>& intensity);
		FILE *create_cel_file (int chipIx);

  /** 
   * @brief Default Getter method for parameters and their documentation.
   * @return map of parameters and their descriptions.
   */
  static std::vector<SelfDoc::Opt> getDefaultDocOptions() 
  { 
    std::vector<SelfDoc::Opt> opts;

	SelfDoc::Opt floor = {"floor", SelfDoc::Opt::Integer, "1", "1", "1", "1024", "probe intensities are floored at this intensity value"};
	  opts.push_back(floor);

	SelfDoc::Opt low = {"low", SelfDoc::Opt::Integer, "20", "20", "2", "1024", "low intensity target value"};
	  opts.push_back(low);

	SelfDoc::Opt high = {"high", SelfDoc::Opt::Integer, "50000", "50000", "2048", "5000000", "high intensity target value"};
	  opts.push_back(high);

	SelfDoc::Opt ceiling = {"ceiling", SelfDoc::Opt::Integer, "1000000", "1000000", "64000", "20000000", "probe intensities are clipped at this intensity value"};
	  opts.push_back(ceiling);

	SelfDoc::Opt low_pct = {"low_pct", SelfDoc::Opt::Double, "0.02", "0.02", "0.001", "0.25", "low target percentile"};
	  opts.push_back(low_pct);

	SelfDoc::Opt high_pct = {"high_pct", SelfDoc::Opt::Double, "0.98", "0.98", "0.75", "0.999", "high target percentile"};
	  opts.push_back(high_pct);

	SelfDoc::Opt gain_lock_pct = {"gain_lock_pct", SelfDoc::Opt::Double, "0.0", "0.0", "0.0", "0.999", "Transition from power series transformation to linear transformation at this intensity percentile."};
	  opts.push_back(gain_lock_pct);

	SelfDoc::Opt cel_out = {"cel_out", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", "true to create a cel file containing Scaled probe intensities This cel file may be used to evaluate the intermediate probe intensities."};
	  opts.push_back(cel_out);


    return opts;
  }

protected:

	ChipLayout *m_ChipLayout_obj;
	bool is_filtered_probe (int iprobe);

	int m_floor;
	int m_low;
	int m_high;
	int m_ceiling;
	double m_gain_lock_pct;

	double m_low_pct;
	double m_high_pct;

	bool m_cel_out;


#if 0
		
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
  //static void learnParameters(const std::vector<float> &data, 
  //                            int binCount, 
  //                            std::vector<vector<float> > &chipBins,
  //                            /*std::vector<int> &probes, */
  //                            std::vector<char> &probeGcVec, 
  //                            bool warnings);

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

  void newChip(std::vector<float> &data);



  /** 
   * Setting the GC counts for all the probes on the array.
   * @param vec - GC count for every probe on array indexed by id.
   */
  void setProbeGcVec(const std::vector<char> &vec);

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

#endif

};

#endif /* _SCIR_H_ */