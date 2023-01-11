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
 * @file   GcCorrection.cpp
 * @author Vincent Bressler
 * @date   Feb 12, 2014
 * 
 * @brief Class for doing an background adjustment on all probes based
 * on the median intensity of probes with similar GC content.
 */


//
#include "chipstream/GcCorrection.h" 
//
#include "chipstream/DataStore.h"
#include "chipstream/QuantMethod.h"
//
#include "stats/stats.h"
#include "util/Err.h"
#include "util/Verbose.h"
#include "util/md5sum.h"
//
#include <vector>
//

#include "AGCC_FileUtilities.h"
#include <iostream>

#ifndef Max
#define Max(a,b) ((a < b) ? b : a)
#endif

#ifndef Min
#define Min(a,b) ((a < b) ? a : b)
#endif      

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif


inline int float_compare( const void *vf1, const void *vf2)
{
	float *f1 = (float *)vf1;
	float *f2 = (float *)vf2;

	if (*f1 < *f2)
		return -1;

	if (*f1 > *f2)
		return 1;

	return 0;
}

inline float PercentileValue (const float *values, const int numValues, const float percentile, float *range=0)
{
	unsigned long valueBefore,valueAfter;
	float percentileIndex,fraction;

	percentileIndex=((float)numValues-1.0F)* percentile;
	valueBefore =(unsigned long)percentileIndex;
	valueAfter = valueBefore + 1;

	// Bounds check on indices
	if (valueAfter > (unsigned long)(numValues-1))
		valueAfter = (unsigned long)(numValues-1);

	if (valueBefore < 0)
		valueBefore = 0;

	fraction = percentileIndex - valueBefore;

	float return_value = values[valueBefore] + fraction * (values[valueAfter] - values[valueBefore]);

	if (range)
		*range = values[valueAfter] - values[valueBefore];

	return return_value;
}

inline float lookup_percentile (float v, const float *sorted_array, int N)
{
	// Not items yet, insert at beginning
	if (N == 0)
		return 0.0F;

	// Greater or equal to the last item, insert at the end
	if (v > sorted_array[N-1])
		return 1.0F;

	// Less than the first item, insert at the beginning
	if (v < sorted_array[0])
		return 0.0F;

	// If the size is less than 5, use a linear search
	int  i;

	if (N < 6)
	{
		for (i = 1; i < N-1; i++)
		{
			if (v <  sorted_array[i])
				return (float)i/(float)N;
		}

		return (float)i/(float)N;
	}

	// Use a binary search
	int istart,iend;

	istart = 0;
	iend = N-1;
	i = N / 2;

	while (1)
	{
		if (v < sorted_array[i])
			iend = i - 1;
		else
			istart = i + 1;

		if (iend - istart < 4)
		{
			for (i = istart; i <= iend; i++)
			{
				if (v < sorted_array[i])
					return (float)i/(float)N;
			}

			return (float)i/(float)N;
		}

		i = (iend - istart) / 2;
		i += istart;
	}
}

using namespace std;
using namespace affx;


void GcCorrection::compute_gcc_arrays (std::vector<float>& data)
{
	clear_gcc_arrays ();

	int invalid_bin_count = 0;

	int bin = 0;

	// Compute the number of probes that have each particular GC Count
	m_probe_count_per_gcc = new int[m_MaxGc+1];
	memset (m_probe_count_per_gcc,0,sizeof(int)*(m_MaxGc+1));
	for (int probeIx = 0; probeIx < data.size (); probeIx++)
	{
		bin = m_BgProbeGc[probeIx];

		if (bin < 0 || bin > m_MaxGc)
		{
			invalid_bin_count++;
			continue;
		}

		m_probe_count_per_gcc[bin]++;
	}

	// Allocate memory to hold the probe intensities in separate arrays, each array containing all of the probes of a particular GC Count
	m_probe_intensities_per_gcc = new float*[m_MaxGc+1];
	memset (m_probe_intensities_per_gcc,0,sizeof(float)*(m_MaxGc+1));
	for (bin = 0; bin <= m_MaxGc; bin++)
	{
		m_probe_intensities_per_gcc[bin] = new float[m_probe_count_per_gcc[bin]];
		memset (m_probe_intensities_per_gcc[bin],0,m_probe_count_per_gcc[bin]*sizeof(float));
	}

	// Assign probe intensities to separate arrays, each array contains all the probes of a particular GC Count

	int *n_per_gcc = new int[m_MaxGc+1];
	memset (n_per_gcc,0,sizeof(int)*(m_MaxGc+1));

	for (int probeIx = 0; probeIx < data.size (); probeIx++)
	{
		bin = m_BgProbeGc[probeIx];

		if (bin < 0 || bin > m_MaxGc)
		{
			invalid_bin_count++;
			continue;
		}

		m_probe_intensities_per_gcc[bin][n_per_gcc[bin]++] = data[probeIx];
	}

	// Sort the probe intensities in separate arrays, each array containing all of the probes of a particular GC Count
	for (bin = 0; bin <= m_MaxGc; bin++)
	{
		qsort (m_probe_intensities_per_gcc[bin],m_probe_count_per_gcc[bin],sizeof(float),float_compare);
	}
}

void GcCorrection::learnParameters(const std::vector<float> &data, int binCount, std::vector<vector<float> > &chipBins,
                           /*std::vector<int> &probes,*/ std::vector<char> &probeGcVec, bool warnings) 
{
	unsigned int probeIx = 0, binIx = 0;
	vector<vector<float> > gcBins(binCount);
	if(probeGcVec.empty()) 
	{
		Err::errAbort("GcCorrection::learnParameters() - Must have some gc control probes set to calculate background.");
	}

	int invalid_bin_count = 0;

	// Fill in the bins with raw data processed by chip streams.
	for(probeIx = 0; probeIx < data.size (); /*probes.size();*/ probeIx++) 
	{
		float intensity = 0;
		unsigned int bin = 0;
//		int probeIndex = probes[probeIx];

		//if(probeGcVec[probeIndex] == (char)NULLPROBEGC)
		//	Err::errAbort("GcBg - Need probe index: " + ToStr(probeIndex) + " GC count for background, but wasn't loaded.");

//		bin = probeGcVec[probeIndex];
		bin = probeGcVec[probeIx];

		if(bin >= binCount || bin < 0) 
		{
			invalid_bin_count++;
			continue;
		}

//		intensity = data[probeIndex];
		intensity = data[probeIx];
		gcBins[bin].push_back(intensity);
	}

	vector<float> binMeds(binCount, 0);
		
	// Estimate the medians for each bin.
	for(binIx = 0; binIx < gcBins.size(); binIx++) 
	{
		// Unlike GcAdjust, we will default to 0 rather than error. Otherwise
		// one would have to have the GC count for all probes on the array
		float med =  0;
		if(gcBins[binIx].empty()) 
		{
			if(warnings)
			Verbose::out(2, "Warning: GC Bin " + ToStr(binIx) + " has no data.");
		}
		else
		{
			med = median_in_place(gcBins[binIx].begin(), gcBins[binIx].end());
		}

		binMeds[binIx] = med;
	}

	chipBins.push_back(binMeds);
}

inline float __log2 (float v)
{
	return log10 (v) * 3.321928F;
}

float GcCorrection::transform(int probeIx, int chipIx, float intensity, int gcBin) 
{
	if(gcBin == (char)NULLPROBEGC) 
	{
		return intensity;
	} 
	else 
	{
		int target_bin = 12;

		if (gcBin == target_bin)
			return intensity;

		float intensity_log2 = __log2(intensity);
		float intensity_bg_log2 = __log2 (m_Bins[chipIx][gcBin]);
		float intensity_bg_target_log2 = __log2(m_Bins[chipIx][target_bin]);

		float adjustment_log2 = intensity_bg_target_log2 - intensity_bg_log2;
		intensity_log2 += adjustment_log2;

		if (intensity_log2 < 1.0F)
			intensity = 2.0F;
		else
			intensity = pow (2,intensity_log2);
	}

	return intensity;
}

#define GENERATE_CEL_FILES

float GcCorrection::transform(int probeIx, int chipIx, float intensity) 
{
	// Look up probe gc count.
	char gc = m_BgProbeGc[probeIx];
	float result = transform(probeIx, chipIx, intensity, gc);
	return result;
}

bool GcCorrection::is_filtered_probe (int iprobe)
{
	if (m_ChipLayout_obj == 0)
		return false;

	if (m_ChipLayout_obj->getQCProbeMap ()[iprobe])
		return true;

	return false;
}

//inline bool is_in_list (int i)
//{
//	if (i == 3934539) return true;
//	if (i ==  538928) return true;
//	if (i == 1469074) return true;
//	if (i ==  619599) return true;
//	if (i == 1495980) return true;
//	if (i == 4896938) return true;
//	if (i == 2056288) return true;
//	if (i ==  818179) return true;
//
//	return false;
//}

void GcCorrection::transform(int chipIx, std::vector<float>& intensity) 
{
	bool bRev2 = false;
	bool bRev4 = true;

	if (m_algo_rev == 4)
	{
		bRev4 = true;
		bRev2 = false;
	}
	else
	if (m_algo_rev == 1)
	{
		bRev4 = false;
		bRev2 = true;
	}

	char cel_file_name_GC[2048];
	memset (cel_file_name_GC,0,sizeof(cel_file_name_GC));

	while (m_cel_out)
	{
		const std::string cel_file_name = m_TransformedIMart->getCelFileNames()[chipIx];


		if (cel_file_name.length ()+8 > sizeof(cel_file_name_GC))
			break;

		const char *c = cel_file_name.c_str ();
		int i=0;
		for (i=0; i<cel_file_name.length (); i++)
			cel_file_name_GC[i] = c[i];

		int len = strlen (cel_file_name_GC);
		if (cel_file_name_GC[len-4] == '.')
			cel_file_name_GC[len-4] = 0;

		if (bRev2)
			strcat (cel_file_name_GC,".GC.cel");
		else if (bRev4)
			strcat (cel_file_name_GC,".GC4.cel");
			
		if (copy_file (cel_file_name.c_str (),cel_file_name_GC) == false)
			memset (cel_file_name_GC,0,sizeof(cel_file_name_GC));

		break;
	}

	if (chipIx == 0)
	{
		m_probe_intensities_per_gcc = 0;
		m_probe_count_per_gcc = 0;
	}

	compute_gcc_arrays (intensity);

	const int gc_count_target = 12;
	int nfeatures_on_the_floor_after = 0;

	for (uint32_t probeIx = 0; probeIx < intensity.size(); probeIx++) 
	{
		if (is_filtered_probe (probeIx) == true)
			continue;
	
		if (bRev2)
			intensity[probeIx] = transform(probeIx, chipIx, intensity[probeIx]);
		else
		if (bRev4)
		{
			//FILE *stream = 0;
			//if (is_in_list (probeIx))
			//{
			//	const char *file_name = "C:\\Users\\Public\\18678166.txt";
			//	bool bExists = access (file_name,00)==0?true:false;
			//	fopen_s (&stream,file_name,"a+");
			//	if (stream && bExists == false)
			//		fprintf (stream,"probeIx\tGCCount\tprobe_count_per_gcc\tpercentile\tprobe_count_per_gcc[12]\treference_intensity\n");
			//}

			int gc_count = m_BgProbeGc[probeIx];
			//if (stream) fprintf (stream,"%ld\t%ld\t",probeIx,gc_count);

			if (gc_count >= 0 && gc_count < m_MaxGc)
			{
				int n = m_probe_count_per_gcc[gc_count];
				//if (stream) fprintf (stream,"%ld\t",n);

				float percentile = lookup_percentile (intensity[probeIx],m_probe_intensities_per_gcc[gc_count],n);
				//if (stream) fprintf (stream,"%f\t",percentile);

				int n_gc_target = m_probe_count_per_gcc[gc_count_target];
				//if (stream) fprintf (stream,"%ld\t",n_gc_target);
				float reference_intensity = PercentileValue (m_probe_intensities_per_gcc[gc_count_target],n_gc_target,percentile);
				//if (stream) fprintf (stream,"%f\n",reference_intensity);

				intensity[probeIx] = reference_intensity;

				if (intensity[probeIx] < 2.0F)
				{
					intensity[probeIx] = 2.0F;
					nfeatures_on_the_floor_after++;
				}
			}

			//if (stream) fclose (stream);
		}
	}

	if (cel_file_name_GC[0])
	{
		char err_msg[1024];
		if (Fusion_UpdateCelFileIntensityValues (cel_file_name_GC,"Default Group",intensity,intensity.size (),err_msg,sizeof(err_msg)) == FALSE)
			unlink (cel_file_name_GC);
	}

	clear_gcc_arrays ();
}

void GcCorrection::newChip(std::vector<float> &data) 
{
	learnParameters(data, m_MaxGc, m_Bins, m_BgProbeGc, m_ChipCount == 0);
}

void GcCorrection::newDataSet(IntensityMart* iMart) 
{
	int dataSetCount = iMart->getCelDataSetCount();

	Verbose::progressBegin(1, "Applying GCCorrection to " + ToStr(dataSetCount) + " cel datasets", dataSetCount, 0, dataSetCount);

	m_TransformedIMart = iMart;
	// if there aren't any chipstream nodes after this, then don't store
	// all of the intensities
	if (m_Streams.empty()) 
	{
		m_TransformedIMart->setStoreAllCelIntensities(false);
	}

	std::vector<float> data;
	for (int d = 0; d < dataSetCount; d++) 
	{
		Verbose::progressStep(1);

		data = iMart->getCelData(d);
		newChip(data);
		m_ChipCount++;
		transform(d, data);
		m_TransformedIMart->setProbeIntensity(d, data);
	}
	Verbose::progressEnd(1, "Done.");

	chipStreamPassNewChip(m_TransformedIMart);
}

void GcCorrection::setControlProbes(std::vector<int> &vec) 
{
	md5sum md5;

	for(int i = 0; i < vec.size(); i++) 
	{
		md5.update_nbo(vec[i]);
	}

	md5.final(m_ProbeMd5Sum);
	setOptValue("subsetmd5", m_ProbeMd5Sum);
	m_BgProbes = vec;
}

void GcCorrection::setControlProbes(std::vector<Probe *> &vec) 
{
	m_BgProbes.resize(vec.size());

	for(int i = 0; i < vec.size(); i++) 
	{
		m_BgProbes[i] = vec[i]->id;
	}

	setControlProbes(m_BgProbes);
}

/** 
 * Setting the GC counts for all the probes on the array.
 * @param vec - GC count for every probe on array indexed by id.
 */
void GcCorrection::setProbeGcVec(const std::vector<char> &vec, ChipLayout *obj/*=0*/) 
{
	m_BgProbeGc.resize(vec.size());
  
	for(int i = 0; i < m_BgProbeGc.size(); i++) 
	{
		m_BgProbeGc[i] = (char) vec[i];
	}

	///////////////////////////////////////////////////////////////////////////////////////
	// Use ChipLayout to exclude control probe sets
	int n_already_not_used = 0;
	int n_control_probes = 0;
	m_ChipLayout_obj = 0;

	if (obj)
	{
		const std::vector<char> qc_map = obj->getQCProbeMap ();

		if (qc_map.size () == m_BgProbeGc.size ())
		{
			m_ChipLayout_obj = obj;

			for (int iprobe=0; iprobe<m_BgProbeGc.size (); iprobe++)
			{
				if (qc_map[iprobe])
				{
					int8_t old_gc_count = m_BgProbeGc[iprobe];
					if (old_gc_count < 0)
						n_already_not_used++;

					n_control_probes++;
				}
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////
}

/** 
 * Setting the GC counts for all the probes on the array.
 * @param vec - GC count for every probe on array indexed by id.
 */
void GcCorrection::setProbeGcVec(const std::vector<unsigned char> &vec) 
{
	if(m_BgProbes.empty()) 
	{
		APT_ERR_ABORT("GcBG::setProbeGcVec() - Must specify background probe ids first.");
	}

	m_BgProbeGc.resize(vec.size());
	for(int i = 0; i < m_BgProbeGc.size(); i++) 
	{
		m_BgProbeGc[i] = (char) vec[i];
	}
}

/**
 * Setup the background probes and background gc indexes by reading
 * them from the blackboard
 */
void GcCorrection::setParameters(PsBoard &board) 
{
	std::vector<int> vec;
	board.getProbeInfo()->getGcControlProbes(vec);
	setControlProbes(vec);
	std::vector<char> gcVec;
	board.getProbeInfo()->getProbeGc(gcVec);
	setProbeGcVec(gcVec);
}
