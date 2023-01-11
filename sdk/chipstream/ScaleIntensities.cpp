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
 * @file   ScaleIntensities.cpp
 * @author Vincent Bressler
 * @date   May 20, 2014
 * 
 * @brief Code to scale cel intensities.
 */

#include "chipstream/ScaleIntensities.h" 
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

inline int roundInt (double v)
{
	return (int)(v+0.5);
}

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
	valueAfter=min(valueAfter,(unsigned long)(numValues-1));
	valueBefore=max(valueBefore,(unsigned long)0);

	fraction = percentileIndex - valueBefore;

	float return_value = values[valueBefore] + fraction * (values[valueAfter] - values[valueBefore]);

	if (range)
		*range = values[valueAfter] - values[valueBefore];

	return return_value;
}



bool examine_intensities (const float *intensities, int ncols, int nrows, double pct_low, double pct_high, double floor, double ceiling, double &v_pct_low, double &v_pct_high, char *emsg, int emsg_alloc_sz)
{
	int N = ncols*nrows;

	float *tmp = new float[N];
	memmove (tmp,intensities,N*sizeof(float));

	qsort ((void*)tmp,N,sizeof(float),float_compare);

	v_pct_low = PercentileValue (tmp,N,(float)pct_low);
	v_pct_high = PercentileValue (tmp,N,(float)pct_high);

	delete [] tmp;

	return true;
}

inline double compute_linear_intensity (double p, double vlow, double v, double target_intensity_low)
{
	double linear_v = pow (p,(v-vlow)) * target_intensity_low;
	return linear_v;
}

double compute_power_for_this_range (double target_range, double vlow, double vhigh, double target_intensity_low)
{
	// Start searching at 2.0
	double p = 2.0;

	double e1 = target_range - compute_linear_intensity (p-.0005, vlow, vhigh, target_intensity_low);
	double e2 = target_range - compute_linear_intensity (p+.0005, vlow, vhigh, target_intensity_low);

	double e = (e2+e1)/2.0;

	double dpde = 0.001 / (e2-e1);

	int iter = 0;

	while (fabs (e) > 0.1)
	{
		iter++;

		p += -e*dpde;

		e1 = target_range - compute_linear_intensity (p-.0005, vlow, vhigh, target_intensity_low);
		e2 = target_range - compute_linear_intensity (p+.0005, vlow, vhigh, target_intensity_low);

		dpde = 0.001 / (e2-e1);

		e = (e2+e1)/2.0;
	}

	return p;
}


double compute_percentile_fast(const float *v, int N, double pct, int approx_to_use)
{
	int n = N;
	int step = 1;
	if (N/approx_to_use >= 2)
	{
		step = N / approx_to_use;

		n = N / step;
	}

	float *tmp = new float[n+1];
	int j = 0;
	for (int i = 0; i < N; i += step)
	{
		tmp[j++] = v[i];

		if (j == n)
			break;
	}
	n = j;

	qsort(tmp, n, sizeof(float), float_compare);

	double v_pct = (double)PercentileValue (tmp,n,(float)pct);
	delete[] tmp;

	return v_pct;
}

inline double compute_slope(double v, double p, double vlow, double target_intensity_low)
{
	double delta = 0.001;

	double v0 = compute_linear_intensity(p, vlow, v - delta, target_intensity_low);
	double v1 = compute_linear_intensity(p, vlow, v + delta, target_intensity_low);

	double slope = (v1 - v0) / (exp(v + delta) - exp(v - delta));

	return slope;
}

bool compute_line_at_transition_point(double p, double vlow, double target_intensity_low, double v_transition, double &slope_after_transition, double &intercept_after_transtion)
{
	double slope = compute_slope(v_transition, p, vlow, target_intensity_low);

	slope_after_transition = slope;

	double x = exp (v_transition);
	double y = compute_linear_intensity(p, vlow, v_transition, target_intensity_low);

	intercept_after_transtion = y - slope*x;

	return true;
}

bool update_intensities (double target_intensity_low, double target_intensity_high, double vlow, double vhigh, double gain_lock_pct, float *intensity_array, int ncols, int nrows, double ceiling, double floor, int &nfloor, int &nceil, char *emsg, int emsg_alloc_sz, const char *histogram_file/*=0*/)
{
	// Find the appropriate power to use for scaling
	double target_range = target_intensity_high-target_intensity_low;
	double p = compute_power_for_this_range (target_range,vlow,vhigh,target_intensity_low);

	double _vmax=-1e20,_vmin=1e20;

	nfloor = 0;
	nceil = 0;

	int N = ncols*nrows;
	
	double v0_re_norm = 0, v1_re_norm = 0;
	double v_transition = 0;
	double slope_after_transition=0;
	double intercept_after_transition = 0;

	bool btransition_computed_successfully = false;

	//double gain_lock_pct = 0.9;
	double re_norm_pct = 0.0;

	if (re_norm_pct)
	{
		v0_re_norm = compute_percentile_fast(intensity_array, N, re_norm_pct, 100000);
		v0_re_norm = pow(2, v0_re_norm);
	}
	
	if (gain_lock_pct)
	{
		v_transition = compute_percentile_fast(intensity_array, N, gain_lock_pct, 100000);
		btransition_computed_successfully = compute_line_at_transition_point(p, vlow, target_intensity_low, v_transition, slope_after_transition, intercept_after_transition);
	}

	for (int idx=0; idx<N; idx++)
	{
		double v = intensity_array[idx];

		double linear_v;

		if (btransition_computed_successfully && v >= v_transition)
			linear_v = exp(v)*slope_after_transition + intercept_after_transition;
		else
			linear_v = compute_linear_intensity (p, vlow, v, target_intensity_low);

		intensity_array[idx] = (float)linear_v;
		if (intensity_array[idx] < floor)
		{
			intensity_array[idx] = (float)floor;
			nfloor++;
		}
		else
		if (intensity_array[idx] > ceiling)
		{
			intensity_array[idx] = (float)ceiling;
			nceil++;
		}

		_vmax = max (linear_v,_vmax);
		_vmin = min (linear_v,_vmin);
	}

	// Do intensity re-normalization
	if (re_norm_pct)
	{
		_vmax = -1e20;
		_vmin = +1e20;
		nfloor = 0;
		nceil = 0;

		v1_re_norm = compute_percentile_fast(intensity_array, N, re_norm_pct, 100000);

		double sf = v0_re_norm / v1_re_norm;

		for (int idx = 0; idx < N; idx++)
		{
			double linear_v = intensity_array[idx];
			linear_v *= sf;

			intensity_array[idx] = (float)linear_v;
			if (intensity_array[idx] < floor)
			{
				intensity_array[idx] = (float)floor;
				nfloor++;
			}
			else
			if (intensity_array[idx] > ceiling)
			{
				intensity_array[idx] = (float)ceiling;
				nceil++;
			}

			_vmax = max(linear_v, _vmax);
			_vmin = min(linear_v, _vmin);
		}
	}

	if (histogram_file)
	{
		int N = ncols*nrows;
		double log_factor = 1.0/log10(p);

		double vmin = log_factor*log10(_vmin);
		double vmax = log_factor*log10(_vmax);
		double range = vmax-vmin;
		const int nbins=200;
		int bins[nbins];
		memset (bins,0,sizeof(bins));

		for (int i=0; i<N; i++)
		{
			int ibin = roundInt (nbins * ((log_factor*log10(intensity_array[i]) - vmin) / range));
			if (ibin < 0) 
				ibin = 0;
			if (ibin >= nbins)
				ibin = nbins-1;

			bins[ibin]++;
		}

		FILE *stream = 0;		
		stream = fopen (histogram_file,"w");
		if (stream)
		{
			fprintf (stream,"p=%lf\n",p);
			fprintf (stream,"log_p(v)(CEL)\tcount\n");

			for (int ibin=0; ibin<nbins; ibin++)
			{
				fprintf (stream,"%lf\t%ld\n",range*((double)ibin/(double)nbins) + vmin,(long)bins[ibin]);
			}
			fclose (stream);
		}
	}

	return true;
}

bool convert_to_log_scale (float *intensities, int ncols, int nrows, char *emsg, int emsg_alloc_sz)
{
	int N = ncols*nrows;
	for (int i=0; i<N; i++)
	{
		double v = intensities[i];
		if (v < 1.0e-20)
		{
            // @todo: rewrite with std::string
			sprintf (emsg,"convert_to_log_scale () invalid intensity %f at index=%ld",v,(long)i);
			return false;
		}

		v = log (v);

		intensities[i] = (float)v;
	}

	return true;
}


bool ReScaleIntensities 
(
	float *intensity_array, 
	int ncols, 
	int nrows,
	
	double pct_low,
	double pct_high,
	
	double target_intensity_low,
	double target_intensity_high,
	
	double floor,
	double ceiling,

	int &nfloor,
	int &nceil,

	double gain_lock_pct,

	char *emsg,
	int emsg_alloc_sz,
	
	const char *histogram_file/*=0*/
)
{
	bool status = convert_to_log_scale (intensity_array, ncols, nrows, emsg, emsg_alloc_sz);
	if (status == false)
		return false;

	double v_pct_low;
	double v_pct_high;

	status = examine_intensities (intensity_array, ncols, nrows, pct_low, pct_high, floor, ceiling, v_pct_low, v_pct_high, emsg, emsg_alloc_sz);

	if (status == false)
		return false;

	status = update_intensities (target_intensity_low, target_intensity_high, v_pct_low, v_pct_high, gain_lock_pct, intensity_array, ncols, nrows, ceiling, floor, nfloor, nceil, emsg, emsg_alloc_sz, histogram_file);
	if (status == false)
		return false;

	return true;
}

FILE *ScaleIntensities::create_cel_file (int chipIx)
{
	const std::string cel_file_name = m_TransformedIMart->getCelFileNames()[chipIx];
	
	char file[2048];
	memset (file,0,sizeof(file));

	if (1)
	{
		if (cel_file_name.length ()+8 < sizeof(file))
		{
			const char *c = cel_file_name.c_str ();
			int i=0;
			for (i=0; i<cel_file_name.length (); i++)
			{
				file[i] = c[i];
			}

			strcat (file,".Scaled.txt");
		}
	}
	else
		sprintf (file,"c:\\Users\\Public\\%ld.txt",(long)chipIx);

	FILE *stream = fopen (file,"w");
	if (stream)
	{
		fprintf (stream,"Idx\tIntensity\n");
	}

	return stream;
}

float *get_float_array (std::vector<float>& v)
{
	float *a = new float[v.size ()];
	for (int i=0; i<v.size (); i++)
		a[i] = v[i];
	return a;
}

void ScaleIntensities::transform(int chipIx, std::vector<float>& intensity_vector) 
{
	double pct_low = m_low_pct; // 0.02;
	double pct_high = m_high_pct; // 0.98;
	double target_intensity_low = m_low; //20.0;
	double target_intensity_high = m_high; //50000.0;
	double floor = m_floor; // 3.0;
	double ceiling = m_ceiling; // 1000000.0;
	double gain_lock_pct = m_gain_lock_pct;

	const char *histogram_file = 0;

	int nfloor = 0;
	int nceil = 0;

	char emsg[2048];

	float *intensity_array = get_float_array (intensity_vector); 
	ReScaleIntensities 
	(
		intensity_array, 
		intensity_vector.size (), 
		1,
	
		pct_low,
		pct_high,
	
		target_intensity_low,
		target_intensity_high,
	
		floor,
		ceiling,

		nfloor,
		nceil,

		gain_lock_pct,

		emsg,
		sizeof(emsg),
	
		histogram_file
	);

	for (uint32_t probeIx = 0; probeIx < intensity_vector.size(); probeIx++) 
	{
		if (is_filtered_probe (probeIx))
			continue;
		intensity_vector[probeIx] = intensity_array[probeIx];
	}

	while (m_cel_out)
	{
		char cel_file_name_SI[2048];
		memset (cel_file_name_SI,0,sizeof(cel_file_name_SI));		
		
		const std::string cel_file_name = m_TransformedIMart->getCelFileNames()[chipIx];

		if (cel_file_name.length ()+8 > sizeof(cel_file_name_SI))
			break;

		const char *c = cel_file_name.c_str ();
		int i=0;
		for (i=0; i<cel_file_name.length (); i++)
			cel_file_name_SI[i] = c[i];

		int len = strlen (cel_file_name_SI);
		if (cel_file_name_SI[len-4] == '.')
			cel_file_name_SI[len-4] = 0;

		strcat (cel_file_name_SI,".SI.cel");

		if (copy_file (cel_file_name.c_str (),cel_file_name_SI) == false)
			break;

		char err_msg[1024];
		if (Fusion_UpdateCelFileIntensityValues (cel_file_name_SI,"Default Group",intensity_vector,intensity_vector.size (),err_msg,sizeof(err_msg)) == 0)
			unlink (cel_file_name_SI);

		break;
	}
}



void ScaleIntensities::newDataSet(IntensityMart* iMart) 
{
	int dataSetCount = iMart->getCelDataSetCount();

	Verbose::progressBegin(1, "Applying ScaleIntensities to " + ToStr(dataSetCount) + " cel datasets", dataSetCount, 0, dataSetCount);

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
		transform(d, data);
		m_TransformedIMart->setProbeIntensity(d, data);
	}
	Verbose::progressEnd(1, "Done.");

	chipStreamPassNewChip(m_TransformedIMart);
}

bool ScaleIntensities::is_filtered_probe (int iprobe)
{
	if (m_ChipLayout_obj == 0)
		return false;

	if (iprobe < 0)
		return false;
	
	if (m_ChipLayout_obj->getQCProbeMap ().size () <= iprobe)
		return false;

	if (m_ChipLayout_obj->getQCProbeMap ()[iprobe])
		return true;

	return false;
}
