////////////////////////////////////////////////////////////////
//
// Copyright (C) 2009 Affymetrix, Inc.
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
 * @file   test-regression-copynumber-familial.cpp
 * @brief  Program for doing regression tests on cyto data.
 */

// This has to be last. Otherwise windows compile fails
#include "calvin_files/utils/src/Calvin.h"
#include "file5/File5.h"
#include "util/Fs.h"
#include "util/FsTestDir.h"
#include "util/LogStream.h"
#include "util/Util.h"
#include "util/Verbose.h"
//
#include <cstring>
#include <string>
//
using namespace std;

class test_regression_copynumber_familial_quick
{
public:
    int numPassed;
    int numFailed;
    std::string testDir; 
    std::set<std::string> setIgnore;
    std::map<std::string, float> mapEpsilon;
    std::vector<std::string> vPaternityFamFile;
    std::vector<std::string> vPaternityFamFile_CytoScanHD;

    test_regression_copynumber_familial_quick()
    {
        numPassed = 0;
        numFailed = 0;

        // Header paramters to ignore
        setIgnore.insert("FileCreationTime");
        setIgnore.insert("FileIdentifier");
        setIgnore.insert("program-version");
        setIgnore.insert("create_date");
        setIgnore.insert("create-date");
        setIgnore.insert("affymetrix-algorithm-param-option-verbose");
        setIgnore.insert("affymetrix-algorithm-param-option-exec-guid");
        setIgnore.insert("affymetrix-algorithm-param-option-program-cvs-id");
        setIgnore.insert("affymetrix-algorithm-param-option-version-to-report");
        setIgnore.insert("affymetrix-algorithm-param-option-command-line");
        setIgnore.insert("affymetrix-algorithm-param-option-mem-usage");
        setIgnore.insert("affymetrix-algorithm-param-option-run-probeset-genotype");
        setIgnore.insert("affymetrix-algorithm-param-option-cels");
        setIgnore.insert("affymetrix-algorithm-param-option-out-dir");
        setIgnore.insert("affymetrix-algorithm-param-state-time-start");
        setIgnore.insert("affymetrix-algorithm-param-state-free-mem-at-start");
        //setIgnore.insert("affymetrix-algorithm-param-option-temp-dir");

        // Data Set columns to ignore (Specify as Group.Set.Column)
//        setIgnore.insert("AlgorithmData.MarkerABSignal.SCAR");

        // Data Set column epsilon override (Specify as Group.Set.Column)
//        mapEpsilon.insert(std::pair<std::string, float>("AlgorithmData.MarkerABSignal.SCAR", (float)0.0005));

        // File tags to check for equivalency
        vPaternityFamFile.push_back("Trio-UPD_UC_U2");
        vPaternityFamFile_CytoScanHD.push_back("Trio-NA32.v3-3_RefSet");
    }

    void execute(const std::string& str) {
      std::string cmd = Fs::convertCommandToUnc(str);
      if (system(cmd.c_str()) != 0) {Err::errAbort("Execution failed for command: " + cmd);}
    }

    void equivalency(const std::string& name, std::vector<string>& vFileTags)
    {
        bool bPassed = true;
        AffxString strFileName1;
        AffxString strFileName2;
        Verbose::out(1, "*");
        for(unsigned int i = 0; (i < vFileTags.size()); i++)
        {

            strFileName1 = "../../../regression-data/data/copynumber/copynumber-familial/" + name + "/" + vFileTags[i] + ".fam";
            strFileName2 = Fs::join(testDir , name , vFileTags[i]) + ".fam";
            if (!Calvin::equivalent(strFileName1, strFileName2, setIgnore, mapEpsilon, 0.0001, 0.9999989, false)) {bPassed = false;}
        }
        Verbose::out(1, "*");
        if (bPassed) {numPassed++;} else {numFailed++;}
    }

    void run()
    {
        paternity();
        paternity_CytoScanHD();
    }

    void paternity()
    {
        string inputPathPrefix = "../../../regression-data/data/copynumber/copynumber-familial/paternity/input/";
        execute("./apt-copynumber-familial \
                    --index-cychp-file " + inputPathPrefix + "UPD_UC-U2a_F_01_NN_20081223.Cytogenetics_Array.REF_MODEL.cychp \
                    --mother-cychp-file " + inputPathPrefix + "UPD_UC-U2b_F_01_NN_20081223.Cytogenetics_Array.REF_MODEL.cychp \
                    --father-cychp-file " + inputPathPrefix + "UPD_UC-U2c_M_01_NN_20081223.Cytogenetics_Array.REF_MODEL.cychp \
                    --allele-frequency-file " + inputPathPrefix + "AlleleFrequencies_Reduced.txt \
                    --analysis paternity.role-validity-threshold=1000.0 \
                    --out-dir " + testDir + "/paternity/output \
                    --familial-file " + vPaternityFamFile[0] + ".fam"
            );
        equivalency("paternity/output", vPaternityFamFile);
    }

    void paternity_CytoScanHD()
    {
        string inputPathPrefix = "../../../regression-data/data/copynumber/copynumber-familial/paternity/input_CytoScanHD/";
        execute("./apt-copynumber-familial \
                    --index-cychp-file " + inputPathPrefix + "NA10839_B12_ReffileR2_CytoScanHD_QY_20101122.cyhd.cychp \
                    --mother-cychp-file " + inputPathPrefix + "NA12006_B7_ReffileR2_CytoScanHD_WC_20101130.cyhd.cychp \
                    --father-cychp-file " + inputPathPrefix + "NA12005_B6_ReffileR2_CytoScanHD_WC_20101130.cyhd.cychp \
                    --allele-frequency-file " + inputPathPrefix + "TrioMarkers.txt \
                    --analysis paternity \
                    --out-dir " + testDir + "/paternity/output_CytoScanHD \
                    --familial-file " + vPaternityFamFile_CytoScanHD[0] + ".fam"
            );
        equivalency("paternity/output_CytoScanHD", vPaternityFamFile_CytoScanHD);
    }
};

int main(int argc, char* argv[]) {

  try {
    FsTestDir testDir;
    testDir.setTestDir("copynumber/copynumber-familial-qt", true);
    
    ofstream logOut;
    string logName;
    logName = Fs::join(testDir.asString(), "test-regression-copynumber-familial.log");
    Verbose::out(1, "Log file: " + logName);
    Fs::mustOpenToWrite(logOut, logName.c_str());
    LogStream log(3, &logOut);

    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);

    Verbose::setLevel(3);
    test_regression_copynumber_familial_quick test;
    test.testDir = testDir.asString();

    if ( !Fs::dirExists(test.testDir + "/paternity/output") ) {
      Fs::mkdirPath(test.testDir + "/paternity/output", false);
    }

    Verbose::out(1, "Execute regression test: test-regression-copynumber-familial");

    test.run();

    Verbose::out(1, "NumPassed: " + ToStr(test.numPassed) + " NumFailed: " + ToStr(test.numFailed) + " for test-regression-copynumber-familial");
    logOut.close();
    return test.numFailed != 0;
  }
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      return 1;
  }
  
  return 1;
}

