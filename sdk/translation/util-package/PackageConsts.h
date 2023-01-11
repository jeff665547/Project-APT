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
// ~/apt2-src/util-package/PackageConsts.h ---
//
// Revision History:
// Created by David Le on 12/17/2014
//

#ifndef __PACKAGECONSTS_H__
#define __PACKAGECONSTS_H__

#include <string>
//------- BATCH FOLDER , AGATHA 1.1 ------------

// --------------- dir names ------------------
#ifndef BAT_AXIOMDATA_DIR
#define BAT_AXIOMDATA_DIR             "AxiomAnalysisSuiteData"
#endif
#ifndef BAT_DATA_DIR
#define BAT_DATA_DIR                  "AnalysisSuiteData"
#endif
#ifndef BAT_SNPOLISHER_DIR
#define BAT_SNPOLISHER_DIR             "SNPolisher"
#endif
#ifndef BAT_QC_DIR
#define BAT_QC_DIR                     "QC"
#endif
#ifndef BAT_OTV_DIR
#define BAT_OTV_DIR                     "OTV_Caller"
#endif

// --------------- file names ------------------

// - - - - - - - - AxiomAnalysisSuiteData folder - - - - - - - - - 
const std::string sBAT_AXIOM_BIN_MASK       = ".bin";
const std::string sBAT_AXIOM_INDEX_TXT_MASK = ".CHP.index.txt";
const std::string sBAT_BLOB_FILE            = "All_genotypes_by_snps.CHP.bin";
const std::string sBAT_BLOB_INDEX_FILE      = "All_genotypes_by_snps.CHP.index.txt";
const std::string sBAT_BATCH_INFO_FILE      = "batch_info.xml";
const std::string sBAT_SAMPLE_INFO_FILE     = "sample_info.txt";
const std::string sBAT_ANALYSIS_SET__FILE   = "AnalysisConfiguration.analysis_settings";
const std::string sBAT_THRESHOLD_SET_FILE   = "AnalysisConfiguration.threshold_settings";

// - - - - - - - - OTV_Caller folder - - - - - - - - - 
const std::string sBAT_OTV_POSTERIORS_FILE  = "OTV.snp-posteriors.txt";

// - - - - - - - - SNPolisher folder - - - - - - - - - 
const std::string sBAT_PS_PERFORMANCE_FILE  = "Ps.performance.txt";
const std::string sBAT_SNP_QC_FILE          = "SNP_QC.txt";

// - - - - - - - - QC folder - - - - - - - - - 
const std::string sBAT_PLATE_QC_SUM_FILE    = "Plate_QCSummary.txt";
const std::string sBAT_PLATE_QC_DETAIL_FILE = "PlateQCDetails.txt";
const std::string sBAT_SAMPLE_QC_SUM_FILE   = "Sample_QCSummary.txt";
const std::string sBAT_GEN_QC_REPORT_FILE   = "GenotypingQC.report.txt";
const std::string sBAT_GEN_QC_RES_FILE      = "Geno-qcResults.txt";
const std::string sBAT_SAMPLE_QC_CEL_FILE   = "sample_QC_cel_files.txt";
const std::string sBAT_SS1_CALLS_FILE       = "AxiomSS1.calls.txt";
const std::string sBAT_SS1_REPORT_FILE      = "AxiomSS1.report.txt";

// - - - - - - - - root folder - - - - - - - - - 
const std::string sBAT_GT_CELL_FILE         = "genotyping_cel_files.txt";
const std::string sBAT_ORIGINAL_CELL_FILE   = "original_cel_files.txt";
const std::string sBAT_GT_REPORT_FILE       = "AxiomGT1.report.txt";
const std::string sBAT_SNP_POSTERIORS_FILE  = "AxiomGT1.snp-posteriors.txt";

// --------------- columns names ------------------
const std::string sBAT_PROBESET_ID = "probeset_id";
const std::string sBAT_CR                   = "CR";
const std::string sBAT_N_MINOR_ALLELE       = "nMinorAllele";
const std::string sBAT_N_AA                 = "n_AA";
const std::string sBAT_N_AB                 = "n_AB";
const std::string sBAT_N_BB                 = "n_BB";
const std::string sBAT_N_NC                 = "n_NC";
const std::string sBAT_HWP_VALUE            = "H.W.p-Value";
const std::string sBAT_MAF                  = "MinorAlleleFrequency";
const std::string sBAT_HW_STATISTICS        = "H.W.statistic";


const int nBIN_FILE_BYTES_PER_SAMPLE = 21;      // Call:byte
                                                //Confidence : float32
                                                // LogRatio : float32
                                                // Strength : float32
                                                // AAlleleSignal : float32
                                                // BAlleleSignal : float32

//------- SUITCASE ZIP , AGATHA 1.0  ------------


#ifndef ZIP_SNPDATAFILE_KEY
#define ZIP_SNPDATAFILE_KEY             "All_genotypes_by_snps.CHP.bin"
#endif
#ifndef ZIP_SNPORIGCALL_KEY
#define ZIP_SNPORIGCALL_KEY             "All_genotypes_by_snps.original.call.bin.index.bin"
#endif
#ifndef ZIP_SNPINDEXFILE_KEY
#define ZIP_SNPINDEXFILE_KEY            "All_genotypes_by_snps.CHP.index.txt"
#endif
#ifndef ZIP_SNPSUMMARYREPORTFILE_KEY
#define ZIP_SNPSUMMARYREPORTFILE_KEY    "Apt_snp_summary_Report.txt"
#endif
#ifndef ZIP_GENOTYPINGREPORTFILE_KEY
#define ZIP_GENOTYPINGREPORTFILE_KEY    "GenotypingReport.txt"
#endif
#ifndef ZIP_PSPERFORMANCEFILE_KEY
#define ZIP_PSPERFORMANCEFILE_KEY       "Ps.performance.txt"
#endif
#ifndef ZIP_SAMPLEINFOFILE_KEY
#define ZIP_SAMPLEINFOFILE_KEY          "SampleInfo.txt"
#endif
#ifndef ZIP_SNPPOSTERIORSFILE_KEY
#define ZIP_SNPPOSTERIORSFILE_KEY       "SnpPosteriorsFile.txt"
#endif
#ifndef ZIP_SUITCASEINFOFILE_KEY
#define ZIP_SUITCASEINFOFILE_KEY        "SuitcaseInfo.xml"
#endif
#ifndef ZIP_CONTENTTYPESFILE_KEY
#define ZIP_CONTENTTYPESFILE_KEY        "[Content_Types].xml"
#endif

#ifndef GUI_INFO
#define GUI_INFO                        "GUI,INFO"
#endif
#ifndef GUI_LOG
#define GUI_LOG                         "GUI,LOG"
#endif
#ifndef GUI_LOGINFO
#define GUI_LOGINFO                     "GUI,INFO,LOG"
#endif

#define OFFSET_CALL                     0   // 1-byte field
#define OFFSET_CONFIDENCE               1   // 4-byte field
#define OFFSET_LOG_RATIO                5   // 4-byte field (1+4)
#define OFFSET_STRENGTH                 9   // 4-byte field (1+4+4)
#define OFFSET_AALLELE_SIGNAL           13  // 4-byte field (1+4+4+4)
#define OFFSET_BALLELE_SIGNAL           17  // 4-byte field (1+4+4+4+4)

#ifndef MAX_LINE
#define MAX_LINE                        65536   // 64KB (2^16)
#endif

#ifndef MAX_BUFFER
#define MAX_BUFFER                      65536   // 64KB (2^16)
#endif

#ifndef __int64
#define __int64                         long long int
#endif

#ifndef __uint64
#define __uint64                        unsigned long long int
#endif

// zlib/minizip constants
#define DEFAULT_GROUP                   0
#define INTENSITY                       0
#define STD_DEV                         1
#define Z_WANT64                        1
#define ZLIB_INTERNAL                   0



#endif
